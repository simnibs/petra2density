# stdlib
import csv
from pathlib import Path
import argparse
import os.path
import sys

# external
from numpy import histogram, argmax, ones, logical_not, polyval
import numpy as np
import nibabel as nib
from nibabel.orientations import aff2axcodes
import matplotlib.pyplot as plt

# simnibs
from simnibs import __version__ as simnibs_version
from samseg import gems
from simnibs.segmentation import charm_main
import simnibs.cli.charm


np.set_printoptions(suppress=True, linewidth=120, precision=6)

MAX_CT_VALUE = 3150 # [hu]
MAX_DENSITY_VALUE = 3147.35469785 # [kg/m3]
DENSITY_WATER = 1000  # [kg/m3]
DENSITY_AIR = 1.275 # [kg/m3]

def cli():
    parser = argparse.ArgumentParser(
                    prog='petra2density',
                    description=(
                        'Takes a T1w and PETRA MRI as input and creates a density map.\n'
                        '\n'
                        'petra2density uses CHARM from SimNIBS to segment and bias-field correct the input images.'
                        '\n'
                        'petra2density accepts arguments for the charm segmentation tool (e.g. --noneck and --forcesform). See charm documentation for details.'
                        ))

    parser.add_argument('subject_id', help="results will be placed in output_folder/subject_id")
    parser.add_argument('t1_path', help="path to the T1w input image. should be .nii or .nii.gz")
    parser.add_argument('petra_path', help="path to the PETRA input image. should be .nii or .nii.gz")
    parser.add_argument('output_folder',
                        help="Folder to store the results. output_folder should already exist. The folder output_folder/subject_id will be created.")

    parser.add_argument('--config_name', default="default",
                        help="Select the config file to use for charm. This should be the name of folder in the config directory. The directory should contain a charm_simnibs.ini file. See the config/default folder for an example.")

    parser.add_argument('--register_to_petra', action='store_true', default=False,
                        help="Registers the T1 image to PETRA before running the rest of the pipeline. This is a good idea when your PETRA image is higher resolution than your T1 image.")

    parser.add_argument('--kplan', action='store_true', default=False,
                        help="Align images to a space compatible with k-Plan (RAS+, origin at LPI).")

    parser.add_argument('--charm_done', action='store_true', default=False,
                        help="use this if charm is already done with the configuration you want to use. if using --charm_done, output_folder/subject_id/m2m_subject_id should already exist.")

    parser.add_argument('--ct_to_density_calibration',
            default="ct_to_density_calibration_cph2025_line_v1.csv",
            choices=[
                "ct_to_density_calibration_cph2025_v1.csv",
                "ct-calibration-low-dose-30-March-2023-v1.csv", # https://github.com/ucl-bug/petra-to-ct/
                ],
            help=("Select the calibration to convert CT to density.\n"
                  "A CSV with points defining a mapping from HU to density. \n"
                  "The file should contain two columns. \n"
                  "First column should be HU values, second column should be density values."))

    parser.add_argument('--norm_petra_to_pct_parameters',
            default="norm_petra_to_pct_parameters_cph2025_v1.csv",
            choices=[
                "norm_petra_to_pct_parameters_cph2025_v1.csv",
                "norm_petra_to_pct_ucl.csv", # https://github.com/ucl-bug/petra-to-ct/
                ],
            help="Select the parameters to convert a normalized petra to pseudo-CT.")

    args, remaining_args = parser.parse_known_args()
    remaining_args.extend([args.subject_id, args.t1_path]) # the t1_path here is not used. it's just here so i can reuse the parseArgument function from charm. the t1_path may be over-written elsewhere in this program
    charm_args = simnibs.cli.charm.parseArguments(remaining_args)
    main(args, charm_args)


def simnibs_version_4_6_or_later():
    v_major, v_minor, _ = [int(v) for v in simnibs_version.split(".")]
    return v_major > 4 or (v_major == 4 and v_minor >= 6)


def charm_settings_file_name():
    return "charm_simnibs_v4-6.ini" if simnibs_version_4_6_or_later() else "charm_simnibs_v4-5.ini"


def main(args, charm_args):
    t1_path = Path(args.t1_path).resolve()
    petra_path = Path(args.petra_path).resolve()
    subject_folder = Path(args.output_folder).resolve() / args.subject_id
    m2m_folder = subject_folder / f"m2m_{args.subject_id}"
    ct_to_density_calibration_path = Path(__file__).parent.resolve() / "maps" / args.ct_to_density_calibration
    norm_petra_to_pct_parameters_path = Path(__file__).parent.resolve() / "maps" / args.norm_petra_to_pct_parameters
    charm_settings = Path(__file__).parent.resolve() / "config" / args.config_name / charm_settings_file_name()

    # read the calibration files. fails here if they are not found
    norm_petra_to_pct_parameters = read_norm_petra_to_pct_parameter_file(norm_petra_to_pct_parameters_path)
    ct_to_density_calibration = read_ct_to_density_calibration_file(ct_to_density_calibration_path)

    if args.charm_done and not m2m_folder.is_dir():
        print(f"m2m_folder is {m2m_folder}")
        print("error. when using --charm_done, the output_folder/subject_id/m2m_subject_id should already exist. It does not. Exiting")
        sys.exit(1)

    if not args.charm_done:
        os.mkdir(subject_folder)
        os.mkdir(m2m_folder)

    if args.kplan and not args.charm_done:
        if args.register_to_petra:
            petra_path = align_image_to_kplan_space(petra_path, subject_folder)
        else:
            t1_path = align_image_to_kplan_space(t1_path, subject_folder)
        
    if args.register_to_petra and not args.charm_done:
        t1_path = register_t1_to_petra(args.subject_id, t1_path, petra_path, subject_folder)

    if not args.charm_done:
        if not simnibs_version_4_6_or_later():
            t1_path = add_small_value_to_image(t1_path, subject_folder)
            petra_path = add_small_value_to_image(petra_path, subject_folder)
        run_segmentation(args.subject_id, t1_path, petra_path, m2m_folder, register_to_petra=args.register_to_petra, use_settings=charm_settings, charm_args=charm_args)

    convert_petra_to_density(m2m_folder, norm_petra_to_pct_parameters, ct_to_density_calibration)


def read_norm_petra_to_pct_parameter_file(path):
    """
    reads the csv file that has the parameters to convert a
    normalized petra to a pseudo-ct

    fails if any parameters are missing and prints out the missing parameters
    """
    params = dict(
            background = None,
            soft_tissue = None,
            bone_offset = None,
            bone_slope = None,
            )

    with open(path) as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if row[0].strip().startswith('#'):
                continue
            params[row[0].strip()] = float(row[1])

    assert not any(v is None for k, v in params.items()), f"must fill all norm_petra_to_pct parameters. missing parameters: {[k for k, v in params.items() if v is None]}"
    return params


def read_ct_to_density_calibration_file(path):
    """
    reads the csv file with the ct to density calibration.

    add a final point with MAX_CT_VALUE and MAX_DENSITY_VALUE. This is to ensure that the maximum threshold is not too low.
    """
    points = np.loadtxt(path, delimiter=",")
    points = np.concatenate((points, [[MAX_CT_VALUE, MAX_DENSITY_VALUE]]))
    hu_values, density_values = points[:, 0], points[:, 1]
    assert np.all(np.diff(hu_values) > 0), "ct to density values must be increasing only in the calibration file."
    assert np.all(np.diff(density_values) > 0), "ct to density values must be increasing only in the calibration file."
    return points


def register_t1_to_petra(subject_id, t1_path, petra_path, output_folder : Path):
    """
    Rigid registration of t1 to petra using charm
    output file is saved as output_folder / t1_reg2petra.nii.gz
    """
    RAS2LPS = np.diag([-1, -1, 1, 1])
    reg = gems.KvlRigidRegistration()
    reg.read_images(str(petra_path), str(t1_path))
    reg.initialize_transform()
    reg.register()
    trans_mat = RAS2LPS@reg.get_transformation_matrix()@RAS2LPS
    t1_output_path = output_folder / "t1_reg2petra.nii.gz"
    reg.write_out_result(str(t1_output_path))
    mat_path = output_folder / 't1_reg2petra_dof6.dat'
    np.savetxt(str(mat_path), trans_mat)
    return t1_output_path


def add_small_value_to_image(image_path, output_folder, small_value=0.01):
    """
    adds a small value to an image and saves it.

    this is sometimes needed for segmenting images with simnibs versions older than 4.6
    """
    image = nib.load(image_path)
    image_data = image.get_fdata()
    min_value = image_data.min()
    if min_value <= 0:
        min_value -= small_value
        print(f"adding {abs(min_value)} to {image_path.name} to ensure all values are positive")
        new_image = nib.Nifti1Pair(image_data + abs(min_value), image.affine, image.header)
        name = image_path.name.replace(".nii", "_positive.nii")
        output_path = output_folder / name
        nib.save(new_image, output_path)
        return output_path
    return image_path


def run_segmentation(subject_id, t1_path, petra_path, m2m_folder, register_to_petra, use_settings, charm_args):
    """
    charm

    """
    charm_main.run(
        str(m2m_folder),
        T1 = str(t1_path),
        T2 = str(petra_path),
        registerT2 = not register_to_petra, # if the t1 has already been registered to the petra then there is no need to register again
        usesettings=str(use_settings),
        initatlas = True,
        segment = True,

        mesh_image = True, # todo: this is on by default for now but not needed for this program
        create_surfaces = True, # todo: this is on by default for now but not needed for this program
        #create_surfaces = charm_args.surfaces,
        #mesh_image = charm_args.mesh,

        noneck = charm_args.noneck,
        init_transform = charm_args.inittransform,
        use_transform = charm_args.usetransform,
        force_qform = charm_args.forceqform,
        force_sform = charm_args.forcesform,
        fs_dir = charm_args.fs_dir,
        options_str = " ".join(sys.argv[1:]),
        debug = charm_args.debug,
    )


def bone_from_label(label):
    """
    get bone label from the m2m_subid/final_tissues.nii.gz label
    """
    return (label == 7) | (label == 8)


def soft_tissue_from_label(label):
    """
    get soft tissue label from the m2m_subid/final_tissues.nii.gz label
    soft tissue is defined as everything that is not air or bone in the image
    """
    background = label == 0
    bone = bone_from_label(label)
    soft_tissue = logical_not(background) & logical_not(bone)
    return soft_tissue


def normalize_petra(petra_data, label, plot=False):
    """
    normalize a petra image to the peak soft-tissue value in the histogram.

    """
    petra = petra_data.squeeze()
    soft_tissue = soft_tissue_from_label(label)
    h = histogram(petra[soft_tissue], bins=100)
    bins = (h[1][1:] + h[1][:-1])/2
    vals = h[0]
    soft_tissue_value = bins[argmax(vals)]
    if plot:
        print(soft_tissue_value)
        plt.plot(bins, vals)
        plt.scatter(soft_tissue_value, max(vals))
        plt.show()
    norm_petra = petra_data / soft_tissue_value
    return norm_petra


def convert_petra_to_density(m2m_folder, norm_petra_to_pct_parameters, ct_to_density_calibration):
    petra_image = nib.load(m2m_folder / "segmentation" / "T2_bias_corrected.nii.gz")
    petra = petra_image.get_fdata()
    
    # save the petra with the name petra
    nib.save(petra_image, m2m_folder / "p2d_petra_bfc.nii.gz")
    
    label_image = nib.load(m2m_folder / "final_tissues.nii.gz")
    label = label_image.get_fdata().squeeze()

    # create and save normalized petra
    norm_petra = normalize_petra(petra, label)
    nib.save(nib.Nifti1Pair(norm_petra, petra_image.affine, petra_image.header), m2m_folder / "p2d_norm_petra.nii.gz")

    # create and save binary bone mask
    bone_mask = bone_from_label(label).astype(np.int16)
    nib.save(nib.Nifti1Pair(bone_mask, petra_image.affine, petra_image.header), m2m_folder / "p2d_bone_mask.nii.gz")

    # create and save pseudo-ct
    pct = petra_to_pct(norm_petra, label, norm_petra_to_pct_parameters)
    nib.save(nib.Nifti1Pair(pct, petra_image.affine, petra_image.header), m2m_folder / "p2d_pct.nii.gz")

    # create and save density
    density = pct_to_density(pct, label, ct_to_density_calibration)
    output_path = m2m_folder / "p2d_density.nii.gz"
    nib.save(nib.Nifti1Pair(density, petra_image.affine, petra_image.header), output_path)
    print(f"Done. Density image is here: {output_path}")
    

def petra_to_pct(norm_petra, label, norm_petra_to_pct_parameters):
    params = norm_petra_to_pct_parameters
    soft_tissue = soft_tissue_from_label(label)
    bone = bone_from_label(label)
    pct = ones(norm_petra.shape) * params["background"]
    pct[soft_tissue] = params["soft_tissue"]
    pct[bone] = params["bone_slope"] * norm_petra[bone] + params["bone_offset"]
    pct[bone & (pct < params["soft_tissue"])] = params["soft_tissue"]
    return pct


def pct_to_density(pct, label, ct_to_density_calibration):
    hu_values, density_values = ct_to_density_calibration[:, 0], ct_to_density_calibration[:, 1]
    density = np.interp(pct[:], hu_values, density_values)
    density = density.reshape(pct.shape)
    bone = bone_from_label(label)
    density[(density < DENSITY_WATER) & bone] = DENSITY_WATER # bone should not have lower density than water
    density[density < DENSITY_AIR] = DENSITY_AIR # nothing should have lower density than air
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
