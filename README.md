# petra2density


`petra2density` converts a T1-weighted MRI and a PETRA MRI of the head into a density image and pCT. Work in progress.


![Example data](/docs/images/example.jpg)


![fitting_hisogram](/docs/images/histogram.png)


See notebooks for examples of data and evaluations steps.


## Requirements

Install SimNIBS. SimNIBS version 4.6 is recommended (4.5 should also work). No other installation steps or additional Python
packages are required.

See the [SimNIBS documentation](https://github.com/simnibs/simnibs) for
installation instructions.


## Installation

Clone this repository and run `petra2density.py` from the SimNIBS Python
environment.

```
git clone https://github.com/simnibs/petra2density.git
```

Use either `simnibs_python` or activate the SimNIBS Conda environment before running the script.


## Input Data

The script expects a T1-weighted NIfTI image and a PETRA NIfTI image of the
same subject.

## Run

Call ``simnibs_python`` or activate the simnibs conda environment.
The ``output_folder`` should already exist. The program creates ``output_folder/SUBJECT_ID/m2m_SUBJECT_ID``.

```
simnibs_python /path/to/petra2density.py \
  SUBJECT_ID \
  /path/to/T1.nii.gz \
  /path/to/PETRA.nii.gz \
  /path/to/output_folder
```

The density image is saved here:

```
/path/to/output_folder/SUBJECT_ID/m2m_SUBJECT_ID/p2d_density.nii.gz
```

## Command-Line Interface

```text
petra2density SUBJECT_ID T1_PATH PETRA_PATH OUTPUT_FOLDER [options] [CHARM options]
```


## Optional arguments


* Use ``--register_to_petra`` to register the T1w image to the PETRA image. This is a good idea if your PETRA image is higher resolution than your T1w image and you want to keep that resolution for your simulations. Default is to register the PETRA image to the T1w image.

* Use ``--kplan`` to align the image to a space compatible with k-plan. Use the kct file in the `resources` folder when running simulation on the density image in k-plan.

* Use ``--ct_to_density_calibration ct_to_density_calibration_cph2025_v1.csv`` with the file name to select the mapping of CT values to density. The file must be in the maps folder. See the ``ct_to_density_calibration_cph2025_v1.csv`` for an example file. The default mapping is ``ct_to_density_calibration_cph2025_v1.csv``.

* Use ``--norm_petra_to_pct_parameters norm_petra_to_pct_parameters_cph2025_v1.csv`` with the file name to select the parameters to use to convert a normalized PETRA image to a pseudo-CT. The file must be in the maps folder. See the ``norm_petra_to_pct_parameters norm_petra_to_pct_parameters_cph2025_v1.csv`` file for an example. The default parameter file is ``norm_petra_to_pct_parameters norm_petra_to_pct_parameters_cph2025_v1.csv``.

* Use `--fill_holes` to fill small soft-tissue holes in bone before creating the pseudo-CT and density images.

## Outputs

Images are saved in `output_folder/SUBJECT_ID/m2m_SUBJECT_ID`.

- `p2d_density.nii.gz`: density image in kg/m3.
- `p2d_pct.nii.gz`: pseudo-CT image in HU.
- `p2d_norm_petra.nii.gz`: Bias-field-corrected PETRA image normalized to the soft-tissue intensity
  peak.
- `p2d_petra_bfc.nii.gz`: Bias-field-corrected PETRA image.
- `p2d_bone_mask.nii.gz`: Binary bone mask.
- `p2d_soft_tissue_bone_label.nii.gz`: Label image with `0` for background,
  `1` for soft tissue, and `2` for bone.
- `final_tissue_with_hole_filling.nii.gz`: Written only when `--fill_holes` is
  used.


## Program steps

1. Read the configuration and parameter files.
2. If run with ``[--kplan]``. Convert the images to a space compatible with k-plan.
3. If run with ``[--register_to_petra]``. Rigidly register the T1 image to the PETRA image using charm.
4. If there are negative values in the images, make the images positive by adding smallest value to the image. This is done to avoid holes in the head mask created by charm. Only done for simnibs versions older than 4.6.
5. Run bias-field correction and head segmentation with charm. Optional parameters for charm can be given to petra2density (e.g. --noneck).
6. Normalize bias-field corrected PETRA image to the peak soft-tissue value in the histogram of soft-tissue voxels.
7. Convert the normalized PETRA image to a pseudo-CT according to the parameters. Can be configured with ``--norm_petra_to_pct_parameters``.
8. Convert the pseudo-CT to density according to the calibration map. Interpolating between points. Can be configured with ``--ct_to_density_calibration``. A max value point of 3150 HU is always added and higher values are thresholded.


## References

Based on work from [petra-to-ct](https://github.com/ucl-bug/petra-to-ct/).
