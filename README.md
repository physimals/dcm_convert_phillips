DCM->Nifti conversion for Phillips data
=======================================

Attempts to convert DICOM files produced by some Phillips scanners using
standard tools (DCM2NIIX, DCMSTACK) did not seem to work, particularly for
4D ASL data.

The underlying reason seems to be difficulty for these tools to determine
the correct ordering of the DICOM slices. This may be caused by non-standard
use of the DICOM format, or may simply reflect ambiguities in the format 
itself.

This script is intended as a fallback conversion method when these tools 
fail. 

**IMPORTANT: If you use this script to convert any of your data, you should
use it to convert all of the data from the same acquisition sequence!**

Usage
-----

    python dcm_convert_phillips.py <dir_name>

Output will be ``<dir_name>.nii``
