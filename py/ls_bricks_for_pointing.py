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
from astropy import wcs

@lru_cache(maxsize=1)
def load_ccd_corners():
    fname = os.path.join(os.environ['DECAM_META'],
                         'cornerCoords_SN-C1-reordered-TAN.fits')

    print('READING ' + fname)

    assert(os.path.exists(fname))

    tab = fits.getdata(fname)

    return tab

@lru_cache(maxsize=1)
def load_brick_wcs_template():
    # would be better to use a pickle file for this
    fname = os.path.join(os.environ['DECAM_META'],
                         'brick_wcs_template-header.fits.gz')

    print('READING ' + fname)

    assert(os.path.exists(fname))

    h = fits.getheader(fname)

    return h

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

def brick_wcs(ra_decam_pointing, dec_decam_pointing):
    # use brick WCS template to make an appropriately
    # centered brick WCS object

    h = load_brick_wcs_template()

    h['CRVAL1'] = ra_decam_pointing
    h['CRVAL2'] = dec_decam_pointing

    w = wcs.WCS(h)

    return w

def decam_wcs():
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

def get_region_name(decam_pointing_dec):
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
