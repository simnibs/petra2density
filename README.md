# petra2density


Convert T1w and PETRA MRI head images to density maps. Work in progress.


![Example data](/docs/images/example.jpg)


# Requirements


* Simnibs 4.5 or higher
* See Simnibs requirements


# Installation


See [the simnibs docs](https://github.com/simnibs/simnibs) for instructions on how to install SimNIBS.


Clone or download the petra2density repository. petra2density can be run from the simnibs environment with no additional dependencies or installation steps required.


# Usage


Call ``simnibs_python`` or activate the simnibs conda environment.


```
simnibs_python /path/to/petra2density.py SUBJECT_ID /path/to/T1.nii.gz /path/to/PETRA.nii.gz /path/to/output_folder
```


The ``output_folder`` should already exist. The program creates ``output_folder/SUBJECT_ID/m2m_SUBJECT_ID``.


The density image will be saved here: ``output_folder/SUBJECT_ID/m2m_SUBJECT_ID/p2d_density.nii.gz``.


## Optional arguments


Use ``--register_to_petra`` to register the T1w image to the PETRA image. This is a good idea if your PETRA image is higher resolution than your T1w image and you want to keep that resolution for your simulations. Default is to register the PETRA image to the T1w image.


Use ``--kplan`` to align the image to a space compatible with k-plan.


Use ``--ct_to_density_calibration ct_to_density_calibration_cph2025_v1.csv`` with the file name to select the mapping of CT values to density. The file must be in the maps folder. See the ``ct_to_density_calibration_cph2025_v1.csv`` for an example file. The default mapping is ``ct_to_density_calibration_cph2025_v1.csv``.


Use ``--norm_petra_to_pct_parameters norm_petra_to_pct_parameters_cph2025_v1.csv`` with the file name to select the parameters to use to convert a normalized PETRA image to a pseudo-CT. The file must be in the maps folder. See the ``norm_petra_to_pct_parameters norm_petra_to_pct_parameters_cph2025_v1.csv`` file for an example. The default parameter file is ``norm_petra_to_pct_parameters norm_petra_to_pct_parameters_cph2025_v1.csv``.


# Program procedure

1. Read the configuration and parameter files.
2. If run with ``[--kplan]``. Convert the images to a space compatible with k-plan.
3. If run with ``[--register_to_petra]``. Rigidly register the T1 image to the PETRA image using charm.
4. If there are negative values in the images, make the images positive by adding smallest value to the image. This is done to avoid holes in the head mask created by charm. Only done for simnibs versions older than 4.6.
5. Run bias-field correction and head segmentation with charm. Optional parameters for charm can be given to petra2density (e.g. --noneck).
6. Normalize bias-field corrected PETRA image to the peak soft-tissue value in the histogram of soft-tissue voxels.
7. Convert the normalized PETRA image to a pseudo-CT according to the parameters. Can be configured with ``--norm_petra_to_pct_parameters``.
8. Convert the pseudo-CT to density according to the calibration map. Interpolating between points. Can be configured with ``--ct_to_density_calibration``. A max value point of 3150 HU is always added and higher values are thresholded.


# References

Based on work from [petra-to-ct](https://github.com/ucl-bug/petra-to-ct/).
