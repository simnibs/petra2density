from numpy import histogram, argmax, ones, logical_not, polyval
import nibabel as nib
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import os.path
import numpy as np


def main():
    parser = argparse.ArgumentParser(
                    prog='petra2density',
                    description='Converts a bias-field-corrected PETRA image to a density map. Requires a CHARM segmentation.')
    parser.add_argument('m2m_folder', help="path to the m2m_folder created by CHARM")
    parser.add_argument('--register_to_petra', action='store_true', default=False)
    parser.add_argument('--ct2density_calibration_file', default=None, help="A CSV with points defining a mapping from HU to density. The file should contain two columns. First column should be HU values, second column should be density values.")
    args = parser.parse_args()
    
    ct2density_calibration_points = np.loadtxt(args.ct2density_calibration_file) if args.ct2density_calibration_file is not None else None    
    
    m2m_folder = Path(args.m2m_folder)

    label_image = nib.load(m2m_folder / "final_tissues.nii.gz")
    label = label_image.get_fdata().squeeze()
    if args.register_to_petra:
        petra_image = nib.load(m2m_folder / "T1.nii.gz")
    else:
        petra_image = nib.load(m2m_folder / "T2_reg.nii.gz")
    petra_data = petra_image.get_fdata()
    
    background = label == 0
    bone = (label == 7) | (label == 8)
    soft_tissue = logical_not(background) & logical_not(bone)
    
    pct = petra_to_pct(bone, soft_tissue, petra_data)
    density = pct_to_density(bone, soft_tissue, pct, ct2density_calibration_points)
    
    density_image = nib.Nifti1Pair(density, petra_image.affine)
    nib.save(density_image, m2m_folder / "density.nii.gz")


def petra_to_pct(bone_label, soft_tissue_label, petra_data, plot=False):
    petra = petra_data.squeeze()
    h = histogram(petra[soft_tissue_label], bins=100)
    bins = (h[1][1:] + h[1][:-1])/2
    vals = h[0]
    soft_tissue_value = bins[argmax(vals)]
    if plot:
        print(soft_tissue_value)
        plt.plot(bins, vals)
        plt.scatter(soft_tissue_value, max(vals))
        plt.show()
    norm_petra = petra_data / soft_tissue_value
    pct = -1000 * ones(norm_petra.shape)
    pct[soft_tissue_label] = 42
    pct[bone_label] = -2929.6 * norm_petra[bone_label] + 3274
    return pct


def pct_to_density(bone_label, soft_tissue_label, ct, ct2density_calibration_points=None):
    rho_water = 1000  # density [kg/m3]
    rho_air = 1.275
    HU_water = 0
    HU_air = -1000
    if ct2density_calibration_points is not None:
        hu_values, density_values = ct2density_calibration_points[:, 0], ct2density_calibration_points[:, 1]
        assert np.all(np.diff(density_values) > 0), "Density values must be increasing only."
        density = np.interp(ct[:], hu_values, density_values)
        density = density.reshape(ct.shape)
        density[(density < rho_water) & bone_label] = rho_water # to avoid having values < water
        density[density < rho_air] = rho_air # otherwise NaN in sensor_data
        return density
    else:
        p = [0.455927942656135, 1.003775211687315e+03]

        density = rho_air * ones(ct.shape)
        density[soft_tissue_label] = rho_water
        density[bone_label] = polyval(p, ct[bone_label]) # HU-rho mapping (HU>0)
        density[(ct < 0) & bone_label] = rho_air + (ct[(ct < 0) & bone_label] - HU_air) * (rho_water-rho_air)/(HU_water-HU_air)  # HU-rho mapping (-1000<HU<0)
        density[(density < rho_water) & bone_label] = rho_water # to avoid having values < water
        density[density < rho_air] = rho_air # otherwise NaN in sensor_data
        return density


if __name__ == "__main__":
    main()
