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
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

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
    # would be better to use a pickle file for this
    fname = os.path.join(os.environ['DECAM_META'],
                         'decam_wcs_template-header.fits.gz')

    print('READING ' + fname)

    assert(os.path.exists(fname))

    h = fits.getheader(fname)

    return h

@lru_cache(maxsize=2)
def load_bricklist(region='south'):
    # use env variable for top-level Legacy Surveys
    # data release directory location
    # call this $LS_ROOT

    # file names:
    # $LS_ROOT/north/survey-bricks-dr9-north.fits.gz
    # $LS_ROOT/south/survey-bricks-dr9-south.fits.gz

    # sort by Dec here for binary search later

    assert(region in ['north', 'south'])

    fname = os.path.join(os.environ['LS_ROOT'], region,
                         'survey-bricks-dr9-' + region + '.fits.gz')

    print('READING ' + fname)
    assert(os.path.exists(fname))

    tab = fits.getdata(fname)

    tab.sort(order='dec')

    return tab

def brick_wcs(racen, deccen):
    # use brick WCS template to make an appropriately
    # centered brick WCS object

    h = load_brick_wcs_template()

    h['CRVAL1'] = racen
    h['CRVAL2'] = deccen

    w = wcs.WCS(h)

    return w

def decam_wcs(ra_decam_pointing, dec_decam_pointing):
    # use DECam tangent plane WCS template to make an appropriately
    # centered DECam tangent plane WCS object
    # does this WCS template run into problems for beyond the pole observing?

    h = load_decam_wcs_template()

    h['CRVAL1'] = ra_decam_pointing
    h['CRVAL2'] = dec_decam_pointing

    w = wcs.WCS(h)

    return w

def get_ccd_corners_radec(ra_decam_pointing, dec_decam_pointing, ccdnum=None):
    # this will use the DECam tangent plane WCS object
    # in combination with the CCD corners table

    tab = load_ccd_corners()
    if ccdnum is not None:
        tab = tab[tab['CCD'] == ccdnum]

    tab = Table(tab)

    w = decam_wcs(ra_decam_pointing, dec_decam_pointing)

    corner_ra, corner_dec = w.all_pix2world(tab['X_TANGENT_PLANE'],
                                            tab['Y_TANGENT_PLANE'], 0)

    tab['corner_ra'] = corner_ra
    tab['corner_dec'] = corner_dec

    return tab

def get_brick_corners_radec(racen, deccen):
    # this will use the brick WCS
    # racen, deccen here are the central coords of the *brick*

    w = brick_wcs(racen, deccen)

    x_corner = [0, 0, 3600, 3600]
    y_corner = [0, 3600, 3600, 0]

    ra, dec = w.all_pix2world(x_corner, y_corner, 0)

    return ra, dec

def trimmed_brick_list_1ccd(ra_corners, dec_corners, region='south'):

    # load full list
    #    trim first by dec range (binary search) then by ang sep

    rect_ccd = SphericalPolygon.from_radec(ra_corners, dec_corners)

    bricks = load_bricklist(region=region)

    dec_max = np.max(dec_corners)
    dec_min = np.min(dec_corners)
    
    inds = np.searchsorted(bricks['DEC'], [dec_min - 0.5, dec_max + 0.5])

    print(len(bricks))

    bricks = bricks[inds[0]:inds[1]]

    if len(bricks) == 0:
        return None

    thresh_deg = 0.5

    sc_ccd = SkyCoord(ra=ra_corners[0]*u.degree, dec=dec_corners[0]*u.degree)

    sc_bricks = SkyCoord(ra=bricks['RA']*u.degree, dec=bricks['DEC']*u.degree)

    ang_sep = sc_ccd.separation(sc_bricks)

    bricks = bricks[ang_sep.deg < thresh_deg]

    if len(bricks) == 0:
        return None

    keep = np.zeros(len(bricks), dtype=bool)
    for i, brick in enumerate(bricks):
        # get ra, dec of brick corners
        ra_brick, dec_brick = get_brick_corners_radec(brick['RA'],
                                                      brick['DEC'])

        # check polygon intersection
        rect_brick = SphericalPolygon.from_radec(ra_brick, dec_brick)
        keep[i] = rect_brick.intersects_poly(rect_ccd)

    bricks = bricks[keep]
    
    return bricks if len(bricks) > 0 else None

def check_poly_overlaps_1ccd(ra_decam_pointing, dec_decam_pointing,
                             ccdnum, region='south', checkplot=True):
    # ccdnum is required?

    # get the ra, dec of CCD corners
    # create a polygon for the CCD
    # get the ra, dec of brick corners (for many bricks)
    #    this involves loading a trimmed version of the brick list
    #    trim first by dec range (binary search) then by ang sep
    
    # create a polygon for each brick
    # check their overlaps

    tab = get_ccd_corners_radec(ra_decam_pointing, dec_decam_pointing,
                                ccdnum=ccdnum)

    bricks = trimmed_brick_list_1ccd(tab['corner_ra'],
                                     tab['corner_dec'])

    if checkplot:
        # guess this will have problems near ra = 0 = 360 ...
        plt.plot(tab['corner_ra'], tab['corner_dec'], c='k')
        #plt.scatter(bricks['RA'], bricks['DEC'], s=20, edgecolor='none')
        for brick in bricks:
            ra_c, dec_c = get_brick_corners_radec(brick['RA'], brick['DEC'])
            plt.plot(ra_c, dec_c)

        plt.xlabel('RA (deg)')
        plt.ylabel('Dec (deg)')

        plt.gca().invert_xaxis()
        plt.show()
                  
    
    return bricks

def get_region_name(decam_pointing_dec):
    # subtleties of Dec boundary still need to be addressed
    # decam_pointing_dec should be in degrees

    # see e.g., https://www.legacysurvey.org/dr9/description/
    # for this threshold value

    dec_thresh = 32.375 # deg
    region = 'north' if decam_pointing_dec > dec_thresh else 'south'

    return region

def get_nearby_bricklist(decam_pointing_ra, decam_pointing_dec, ccdnum=1):
    # decam_pointing_ra, decam_pointing_dec should be in degrees
    # figure out relevant sky region (north or south)
    # probably some subtleties right at the boundary, worry about those later..
    # load the full (Dec-sorted) bricklist for relevant sky region
    # cut down to relevant Dec range, with some padding
    # angular separation computation and cut

    # check polygon overlaps
    # return list of overlapping bricks
    
    pass
