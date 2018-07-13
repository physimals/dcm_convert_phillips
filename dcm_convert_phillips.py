"""
Moderately evil script to convert DICOMS from Phillips scanner to a NIFTI file.

DCMSTACK and dcm2niix do not seem to work because of the ordering and use of non-standard
DICOM tags to describe the sequence. In general if either of these tools is unable
to order the 2D DICOM images correctly, 4D datasets will come out completely wrong
and 3D data sets may appear to be correct but actually have incorrect orientation
information.

This script was designed for ASL data sets with label/control images and possibly
multiple PLD/TI images. However it should also work for 3D images from the same
scanner, e.g. structural images.

IMPORTANT: If you need to use this script to convert any of your DICOMS, you should
also use it to convert the other DICOMS from the same sequence. Otherwise you
may end up with multiple NIFTI files which do not align because of the problems
with determining orientation when the ordering of the slices is non-standard.
"""
import sys
import os
import glob

import numpy as np

import nibabel as nib

# Older versions of pydicom imported under a different name
try:
    import dicom
except ImportError:
    import pydicom as dicom

# Vendor-specific DICOM fields we need to use for ASL ordering
DCM_LBL_CTL = (0x2005, 0x1429)
DCM_PHASE_NUM = (0x2001, 0x1008)

def convert_dicoms(dicom_dir, out_fname=None, log=sys.stdout):
    """
    Convert DICOMS in a directory to a NIFTI file

    :arg dicom_dir: Path the directory containing DICOM files
    :arg out_fname: Output NIFTI filename
    """
    log.write("\nReading DICOMS in %s...   0%%" % dicom_dir)
    dcm_files = glob.glob('%s/*' % dicom_dir)

    # All the DICOM images. 
    # 
    # This is a dictionary mapping values of the label/control field
    # to a dictionary mapping values of the phase number field to a sequence
    # of DICOM objects. This rather convoluted structure enables us to 
    # easily group label/control images and images with matching phase together
    # which would otherwise have identical positioning information to other
    # images. This is only relevant for ASL data sets, for other data sets
    # the label/control attribute is blank and the phase is zero.
    dicom_images = {}

    # List of any files found which could not be parsed as DICOM images
    non_dcm_files = []

    # All phases numbers and slice positions found. Use a set so only
    # distinct values are stored
    phases, slice_positions, lbl_ctl_values = set(), set(), set()

    for idx, dcm_path in enumerate(dcm_files):
        try:
            dcm = dicom.read_file(dcm_path)
        except:
            non_dcm_files.append(dcm_path)
            continue
        
        # Get the phase number and label/control label if present
        is_lbl = DCM_LBL_CTL in dcm and dcm[DCM_LBL_CTL].value == "LABEL"
        lbl_ctl_values.add(is_lbl)

        if DCM_PHASE_NUM in dcm:
            phase = dcm[DCM_PHASE_NUM].value
        else:
            phase = 0
        phases.add(phase)

        # Store the dicom object in the dictionary structure
        dicom_images[is_lbl] = dicom_images.get(is_lbl, {})
        dicom_images[is_lbl][phase] = dicom_images[is_lbl].get(phase, [])
        dicom_images[is_lbl][phase].append(dcm)

        # Record the slice position in our set of distinct positions
        slice_positions.add(float(dcm.SliceLocation))

        # Report progress
        percent = 100*float(idx+1) / len(dcm_files)
        log.write("\b\b\b\b%3i%%" % int(percent))
        log.flush()

    log.write("\b\b\b\bDONE\n\n")

    log.write("Slice locations: %s\n" % ", ".join(["%.1f" % s for s in sorted(slice_positions)]))
    log.write("Phases: %s\n" % ", ".join([str(p) for p in sorted(phases)]))
    log.write("Ignored (non-DICOM) files: %i\n" % len(non_dcm_files))
    for f in non_dcm_files: log.write("   - %s\n" % os.path.basename(f))

    # In this section we update the InstanceNumber attribute on each DICOM
    # object to a unique index which orders the images correctly. 
    #
    # The overall (outer) ordering is by phase number, so in an ASL data
    # set all the images for a particular TI/PLD will be together.
    #
    # The next level of ordering is by repeats, i.e. images having the 
    # phase are ordered using the existing InstanceNumber attribute on 
    # the DICOM object (the same one we will later update)
    #
    # Next, images with the same phase number and repeat index are ordered
    # by label/control value, label images first.
    # 
    # Finally, we are down to the slices making up a single 3D image.
    # These are ordered by the slice position 
    #
    # Each DICOM image is added to a simple list - now we have set a unique
    # InstanceNumber on each we will be able to sort them into the correct order
    log.write("\nDetermining ordering of DICOM slices...")
    log.flush()

    nslices = len(slice_positions)
    dicom_images_indexed = []

    for is_lbl in (True, False):
        if is_lbl not in dicom_images: 
            # We may not have any images flagged as labels
            continue
        for phidx, phase in enumerate(sorted(phases)):
            # How many repeats - i.e. distinct images with the same phase  
            nrepeats = len(dicom_images[is_lbl][phase])/nslices

            for sidx, slice_pos in enumerate(sorted(slice_positions)):   

                # Get all the images with a given slice position
                vols = [img for img in dicom_images[is_lbl][phase] if img.SliceLocation == slice_pos]
                for ridx, vol in enumerate(sorted(vols, key=lambda x: x.InstanceNumber)):
                    # Work out the ordering index, as described above and update the InstanceNumber
                    # attribute on the DICOM to match
                    inum = sidx + nslices*int(is_lbl) + nslices*2*ridx + nslices*2*nrepeats*phidx
                    vol.InstanceNumber = inum
                    dicom_images_indexed.append(vol)

    log.write("DONE\n")

    # We have now figured out what order the slices should be put in to make a 4D volume 
    # and stored it in the InstanceNumber DICOM attribute. Now we create a 4D Numpy array
    # and insert the slice data into it in this order
    log.write("Creating ordered 4D volume...")
    log.flush()

    sorted_vols = sorted(dicom_images_indexed, key=lambda x: x.InstanceNumber)
    dcm_first, dcm_last = sorted_vols[0], sorted_vols[-1]
    shape_2d = dcm_first.pixel_array.shape 
    data = np.zeros([shape_2d[0], shape_2d[1], len(slice_positions), int(len(dicom_images_indexed)/len(slice_positions))])
    sidx, vidx = 0, 0
    for idx, dcm in enumerate(sorted_vols):
        data[:, :, sidx, vidx] = np.squeeze(dcm.pixel_array)
        sidx += 1
        if sidx == len(slice_positions):
            sidx = 0
            vidx += 1

    log.write("DONE\n")

    # Figure out the DCM affine
    # Using the method documented at: http://nipy.org/nibabel/dicom/dicom_orientation.html
    # Note that this has to be done after ordering as it depends on the values in the 'first' and 'last' DCM
    log.write("Determining DICOM orientation transformation...")
    log.flush()

    dcm_affine = np.identity(4)
    dcm_affine[:3, 0] = np.array(dcm_first.ImageOrientationPatient[3:]) * dcm_first.PixelSpacing[0]
    dcm_affine[:3, 1] = np.array(dcm_first.ImageOrientationPatient[:3]) * dcm_first.PixelSpacing[1]
    dcm_affine[:3, 3] = np.array(dcm_first.ImagePositionPatient)
    dcm_affine[:3, 2] = (np.array(dcm_last.ImagePositionPatient) - dcm_affine[:3, 3])/(len(slice_positions)-1)

    # The DICOM coordinate system is left-posterior-superior rather 
    # than right-anterior-superior for NIFTI, so we need to invert the first two axes
    dcm_affine[0, :] = -dcm_affine[0, :]
    dcm_affine[1, :] = -dcm_affine[1, :]
    log.write("DONE\n")
    log.write("%s\n" % dcm_affine)

    # Finally, create a NIFTI file from the 4D Numpy array and the affine we have determined.
    log.write("Creating NIFTI...")
    log.flush()

    nii_out = nib.Nifti1Image(data, dcm_affine)
    nii_out.update_header()
    if out_fname is None:
        if dcm_first.SeriesDescription:
            out_fname = "%s.nii" % dcm_first.SeriesDescription
        else:
            log.write("WARNING: Using default output filename\n")
            out_fname = "dcm_to_nifti.nii"
    nii_out.to_filename(out_fname)

    log.write("DONE\n")
    log.write("Wrote %s.\n" % out_fname)

def main():
    """
    Entry point for command line usage
    """
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: dcmconvert_phillips <dicom dir>\n")
        sys.exit(1)

    fpath = sys.argv[1]
    convert_dicoms(fpath, log=sys.stdout)

if __name__ == "__main__":
    main()
