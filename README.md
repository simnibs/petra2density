# petra2density
Convert PETRA head images to density maps. Work in progress, things may change.

![Example data](/docs/images/example.jpg)

# Requirements

* Linux
* Python
* Simnibs

See [the simnibs docs](https://github.com/simnibs/simnibs) for installation instructions.

# Usage

Activate the simnibs conda environment.


```
python petra2density.py SUBJECT_ID /path/to/T1.nii.gz /path/to/PETRA.nii.gz /path/to/output_folder [--register_to_petra] [--ct2density_calibration_file /path/to/calibration_file.csv] [--kplan]
```


Output_folder should already exist. Creates ``output_folder/SUBJECT_ID/m2m_SUBJECT_ID``.


The density image is here: ``output_folder/SUBJECT_ID/m2m_SUBJECT_IDdensity.nii.gz``.


Use ``--register_to_petra`` to register to the PETRA image. Default is to register to the T1 image.


Use ``--kplan`` to align the image to a space compatible with kplan.

# References

Based on work from [petra-to-ct](https://github.com/ucl-bug/petra-to-ct/).
