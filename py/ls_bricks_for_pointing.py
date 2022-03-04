import pandas as pd
from astropy.table import Table
import numpy as np
from spherical_geometry.polygon import SphericalPolygon
import time
import matplotlib.pyplot as plt
import fitsio
import astropy.io.fits as fits
import os
from functools import lru_cache

# lru_cache here eventually
@lru_cache(maxsize=1)
def load_ccd_corners():
    pass

@lru_cache(maxsize=1)
def load_brick_wcs_template():
    # env variable called DECAM_META

    fname = os.path.join(os.environ['DECAM_META'], 'brick_wcs_template.fits')

    print('READING ' + fname)

    assert(os.path.exists(fname))

    result = fits.getdata(fname)

    return result

@lru_cache(maxsize=1)
def load_decam_wcs_template():
    # env variable called DECAM_META
    
    pass

@lru_cache(maxsize=2)
def load_bricklist(region='south'):
    # use env variable for top-level Legacy Surveys
    # data release directory location
    # call this $LS_ROOT

    # file names that I want:
    # $LS_ROOT/north/survey-bricks-dr9-north.fits.gz
    # $LS_ROOT/south/survey-bricks-dr9-south.fits.gz

    # sort by Dec here for binary search later

    # https://docs.astropy.org/en/stable/table/modify_table.html
    # https://numpy.org/doc/stable/reference/generated/numpy.recarray.sort.html
    
    pass

def brick_wcs_instance():
    # use brick WCS template to make an appropriately
    # centered brick WCS object
    pass

def decam_wcs_instance():
    # use DECam tangent plane WCS template to make an appropriately
    # centered DECam tangent plane WCS object
    pass

def get_ccd_corners_radec():
    # this will use the DECam tangent plane WCS object
    # in combination with the CCD corners table
    pass

def get_brick_corners_radec():
    # this will use the brick WCS template
    pass

def check_poly_overlaps():
    pass

def get_sky_region(decam_pointing_dec):
    # subtleties of Dec boundary still need to be addressed
    # decam_pointing_dec should be in degrees

    # see e.g., https://www.legacysurvey.org/dr9/description/
    # for this threshold value

    dec_thresh = 32.375 # deg
    region = 'north' if decam_pointing_dec > dec_thresh else 'south'

    return region

def get_nearby_bricklist(decam_pointing_ra, decam_pointing_dec):
    # decam_pointing_ra, decam_pointing_dec should be in degrees
    # figure out relevant sky region (north or south)
    # probably some subtleties right at the boundary, worry about those later..
    # load the full (Dec-sorted) bricklist for relevant sky region
    # cut down to relevant Dec range, with some padding
    # angular separation computation and cut

    # check polygon overlaps
    # return list of overlapping bricks
    
    pass
