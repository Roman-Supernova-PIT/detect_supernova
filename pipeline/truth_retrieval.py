import os
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from coord_projection import xy_in_image, radec_in_image, two_direction_skymatch

MATCH_RADIUS = 0.4
IMAGE_WIDTH =4088
IMAGE_HEIGHT = 4088
OFFSET = 0

def merge_science_and_template_truth(science_truth, template_truth, science_wcs, template_wcs,
                                     match_radius=None, image_width=None, image_height=None, offset=None):

    match_radius = match_radius or MATCH_RADIUS
    image_width = image_width or IMAGE_WIDTH
    image_height = image_height or IMAGE_HEIGHT
    offset = offset or OFFSET

    science_in_science = xy_in_image(science_truth.x, science_truth.y, width=image_width, height=image_height, offset=offset)
    template_in_template = xy_in_image(template_truth.x, template_truth.y, width=image_width, height=image_height, offset=offset)

    science_truth = science_truth[science_in_science].copy().reset_index(drop=True)
    template_truth = template_truth[template_in_template].copy().reset_index(drop=True)

    science_in_template = radec_in_image(science_truth.ra, science_truth.dec, wcs=template_wcs, width=image_width, height=image_height, offset=offset)
    template_in_science = radec_in_image(template_truth.ra, template_truth.dec, wcs=science_wcs, width=image_width, height=image_height, offset=offset)
    
    science_truth = science_truth[science_in_template].copy().reset_index(drop=True)
    template_truth = template_truth[template_in_science].copy().reset_index(drop=True)

    science_skycoord = SkyCoord(science_truth.ra, science_truth.dec, frame='icrs', unit='deg')
    template_skycoord = SkyCoord(template_truth.ra, template_truth.dec, frame='icrs', unit='deg')


    matched_status, matched_id = two_direction_skymatch(template_skycoord, science_skycoord, radius=match_radius * u.arcsec)
    template_truth_matched = template_truth[matched_status].copy().reset_index(drop=True)
    template_truth_unmatched = template_truth[~matched_status].copy().reset_index(drop=True)
    template_truth_unmatched['flux'] = template_truth_unmatched['flux'] * -1

    science_truth_new = science_truth.copy()
    science_truth_new.loc[matched_id[matched_status], 'flux'] = science_truth_new.loc[matched_id[matched_status], 'flux'] - template_truth_matched.flux.to_numpy()
    
    merged_truth = pd.concat([science_truth_new, template_truth_unmatched]).reset_index(drop=True)
    merged_truth = merged_truth.drop(['x', 'y', 'mag'], axis=1)
    
    return merged_truth