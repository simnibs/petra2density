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
./petra2density.sh [-p] SUBJECT_ID /path/to/T1.nii.gz /path/to/PETRA.nii.gz
```

Creates ``m2m_SUBJECT_ID`` in the current directory with ``density.nii.gz``.

Use ``-p`` to register to the PETRA image. Default is to register to the T1 image.

# References

Based on work from [petra-to-ct](https://github.com/ucl-bug/petra-to-ct/).
