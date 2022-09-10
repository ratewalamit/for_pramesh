#apply hsc masks
from scipy.spatial import KDTree
from astropy.table import Column,Table
import numpy as np
from astropy.io import fits
import pandas as pd
import healpy,sys
from scipy.spatial import cKDTree

input_file=sys.argv[1]


def apply_hsc_masks(fname):
    mask = healpy.fitsfunc.read_map("/mnt/home/faculty/csurhud/github/weaklens_pipeline/DataStore/S16A_v2.0/S16A_FDFC_with_starmask.fits", dtype=int)
    data = fits.open(fname)[1].data
    ra = data["ra"]
    dec = data["dec"]

    theta = (90-dec)*np.pi/180.
    phi = ra*np.pi/180.
    idx = mask[healpy.pixelfunc.ang2pix(4096, theta, phi)] == 1
    data=data[idx]

    table_save=Table(data)
    applied_masks=f"{fname[:-5]}_in_hsc_fov.fits"
    table_save.write(applied_masks,overwrite=True)
    print(f"Masks applied and file saved with name \"{fname[:-5]}_in_hsc_fov.fits\"")
    return (applied_masks)

if __name__ == "__main__":
    if input_file=="":
        print("Enter file to apply masks")
        exit()
    output_file=apply_hsc_masks(input_file)
    print(output_file)


