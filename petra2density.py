from numpy import histogram, argmax, ones, logical_not, polyval
import nibabel as nib
from nibabel.orientations import aff2axcodes
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import os.path
import numpy as np
import sys
import subprocess

np.set_printoptions(suppress=True, linewidth=120, precision=6)

def cli():
    parser = argparse.ArgumentParser(
                    prog='petra2density',
                    description='Runs a head and brain segmentation with CHARM and converts a bias-field-corrected PETRA image to a density map.')

    parser.add_argument('subject_id', help="results will be placed in output_folder/subject_id")
    parser.add_argument('t1_path', help="path to the T1w input image. should be .nii or .nii.gz")
    parser.add_argument('petra_path', help="path to the PETRA input image. should be .nii or .nii.gz")
    parser.add_argument('output_folder',
                        help="Folder to store the results. output_folder should already exist. The folder output_folder/subject_id will be created.")

    parser.add_argument('--register_to_petra', action='store_true', default=False,
                        help="Swaps the input arguments to CHARM so that T1 images is registered to PETRA")

    parser.add_argument('--ct2density_calibration_file', default=None,
                        help=("A CSV with points defining a mapping from HU to density. "
                              "The file should contain two columns. "
                              "First column should be HU values, second column should be density values."))

    parser.add_argument('--kplan', action='store_true', default=False,
                        help="Align images to a space compatible with k-Plan (RAS+, origin at LPI).")

    parser.add_argument('--charm_done', action='store_true', default=False,
                        help="use this if charm is already done with the configuration you want to use. if using --charm_done, output_folder/subject_id/m2m_subject_id should already exist.")

    args = parser.parse_args()

    t1_path = Path(args.t1_path).resolve()
    petra_path = Path(args.petra_path).resolve()
    subject_folder = Path(args.output_folder).resolve() / args.subject_id
    m2m_folder = subject_folder / f"m2m_{args.subject_id}"
    charm_settings = Path().resolve() / "charm.ini"
    petra_is_called_t1_in_m2m_folder = args.register_to_petra

    if args.charm_done and not m2m_folder.is_dir():
        print(f"m2m_folder is {m2m_folder}")
        print("error. when using --charm_done, output_folder/subject_id/m2m_subject_id should already exist. this is the . It does not. Exiting")
        sys.exit(1)

    if not args.charm_done:
        os.mkdir(subject_folder)

    if args.kplan:
        if args.register_to_petra:
            petra_path = align_image_to_kplan_space(petra_path, subject_folder)
        else:
            t1_path = align_image_to_kplan_space(t1_path, subject_folder)

    if not args.charm_done:
        run_segmentation(args.subject_id, t1_path, petra_path, subject_folder, register_to_petra=args.register_to_petra, use_settings=charm_settings)
    convert_petra_to_density(m2m_folder, args.ct2density_calibration_file, petra_is_called_t1_in_m2m_folder)


def run_segmentation(subject_id, t1_path, petra_path, output_folder, register_to_petra=False, use_settings=None):
    first_image, second_image = t1_path, petra_path
    if register_to_petra:
        first_image, second_image = petra_path, t1_path

    subprocess.run(["charm", subject_id, first_image, second_image, "--usesettings", use_settings, "--forceqform"], cwd=output_folder) 


def convert_petra_to_density(m2m_folder, ct2density_calibration_file=None, petra_is_called_t1=False):
    if ct2density_calibration_file is not None and ct2density_calibration_file != "none":
        print("Using ct->density calibration file: {ct2density_calibration_file}")
        ct2density_calibration_points = np.loadtxt(ct2density_calibration_file)
    else:
        ct2density_calibration_points = None
    
    label_image = nib.load(m2m_folder / "final_tissues.nii.gz")
    label = label_image.get_fdata().squeeze()
    if petra_is_called_t1:
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
    output_path = m2m_folder / "density.nii.gz"
    nib.save(density_image, output_path)
    print(f"Done. Density image is here: {output_path}")
    

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
        density[(density < rho_water) & bone_label] = rho_water # to avoid having values less than water
        density[density < rho_air] = rho_air # no density values should be less than air
        return density
    else:
        p = [0.455927942656135, 1.003775211687315e+03]

        density = rho_air * ones(ct.shape)
        density[soft_tissue_label] = rho_water
        density[bone_label] = polyval(p, ct[bone_label]) # HU-rho mapping (HU>0)
        density[(ct < 0) & bone_label] = rho_air + (ct[(ct < 0) & bone_label] - HU_air) * (rho_water-rho_air)/(HU_water-HU_air)  # HU-rho mapping (-1000<HU<0)
        density[(density < rho_water) & bone_label] = rho_water # to avoid having values less than water
        density[density < rho_air] = rho_air # no density values should be less than air
        return density
 
 
def describe(img, title=""):
    hdr = img.header
    aff = img.affine
    ax = aff2axcodes(aff)
    zooms = hdr.get_zooms()[:3]
    print(f"\n=== {title} ===")
    print(f"shape: {img.shape}")
    print(f"zooms: {zooms}")
    print(f"axcodes: {ax}  (expected: RAS)")
    print("affine (voxel -> world):\n", aff)
 

def corners_world(aff, shape):
    """Return world coords of the 8 volume corners (min/max indices along each axis)."""
    I = np.array([[0,0,0,1],
                  [shape[0]-1, 0, 0, 1],
                  [0, shape[1]-1, 0, 1],
                  [0, 0, shape[2]-1, 1],
                  [shape[0]-1, shape[1]-1, 0, 1],
                  [shape[0]-1, 0, shape[2]-1, 1],
                  [0, shape[1]-1, shape[2]-1, 1],
                  [shape[0]-1, shape[1]-1, shape[2]-1, 1]], dtype=float)
    W = (aff @ I.T).T[:, :3]
    return W
 

def align_image_to_kplan_space(input_path, output_folder):
    """
    Align a T1-weighted NIfTI to k-Plan / dispatch world coordinates:
    - Reorient to RAS+ (axes increasing Right, Anterior, Superior)
    - Shift affine so the *true* LPI corner (min over all corners) is at world origin (0,0,0)
    - Write a new NIfTI with updated sform/qform
    - Print a verification report
    """
    if not os.path.isfile(input_path):
        print(f"ERROR: File not found: {input_path}", file=sys.stderr)
        sys.exit(1)
 
    output_name = os.path.basename(input_path).replace(".nii.gz","").replace(".nii","") + "_kplan.nii.gz"
    out_path = output_folder / output_name
 
    print("\nLoading image...")
    img = nib.load(input_path)
    describe(img, "Original image")
 
    # Step 1: Reorient to RAS+ using NiBabel's canonical function
    print("\nReorienting to RAS+ (if needed)...")
    img_ras = nib.as_closest_canonical(img, enforce_diag=False)
    describe(img_ras, "After canonical (RAS+)")
 
    # Step 2: Shift origin so the TRUE LPI corner (min over all corners) maps to (0,0,0)
    aff = img_ras.affine.copy()
    W = corners_world(aff, img_ras.shape)
    mins = W.min(axis=0)
    print(f"\nCorner-wise mins BEFORE shift (L,P,I): {mins}")
    aff[:3, 3] -= mins  # translate so that min over x,y,z is 0
 
    # Step 3: Create new image with updated affine; set both sform and qform
    hdr = img_ras.header.copy()
    sform_code = 1  # NIFTI_XFORM_SCANNER_ANAT (use 2 for ALIGNED_ANAT if preferred)
    hdr.set_sform(aff, code=sform_code)
    try:
        hdr.set_qform(aff, code=sform_code)
    except Exception as e:
        print(f"Warning: could not set qform from affine ({e}). Proceeding with sform only.")
 
    out_img = nib.Nifti1Image(img_ras.get_fdata(dtype=np.float32), aff, header=hdr)
    nib.save(out_img, out_path)
    print(f"\nSaved: {out_path}")
 
    # Step 4: Verify
    print("\nVerifying output...")
    chk = nib.load(out_path)
    describe(chk, "Output (should be RAS+)")
 
    zero_world_after = (chk.affine @ np.array([0,0,0,1.0]))[:3]
    print(f"World coords of voxel (0,0,0) AFTER shift:  {zero_world_after} (should be ~[0,0,0])")
 
    W_after = corners_world(chk.affine, chk.shape)
    mins_after = W_after.min(axis=0)
    maxs_after = W_after.max(axis=0)
    print("\nWorld-space bounding box (after alignment):")
    print("  min(L,P,I) ~", mins_after)
    print("  max(R,A,S) ~", maxs_after)
 
    # Simple sanity checks
    ax = aff2axcodes(chk.affine)
    ok_axes = (ax == ('R','A','S'))  # boolean, not iterable
    ok_origin = np.allclose(mins_after, np.zeros(3), atol=1e-4)  # all mins ~ 0
    print(f"\nChecks: RAS axes: {ok_axes}; LPI corner at (0,0,0): {ok_origin}")
    if not ok_axes:
        print("WARNING: Axes are not RAS+. Something unusual about the affine/orientation.")
    if not ok_origin:
        print("WARNING: Corner mins are not exactly zero. Review affine or tolerance.")
 
    print("\nDone.")

    return out_path
 

if __name__ == "__main__":
    cli()
