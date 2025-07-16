import os

import numpy as np
import pandas as pd

from detect_supernova.make_openuniverse_subtraction_pairs import get_image_info_for_ra_dec, get_templates_for_points
from detect_supernova.util import get_center_and_corners


def test_query():
    ra, dec = 8.2, -43.0
    images = get_image_info_for_ra_dec(ra, dec, band=None)

    assert len(images) == 1110


def test_get_templates_for_points():
    band = "R062"

    corners = [
        (8.2, -43.1),
        (8.2, -42.9),
        (8.25, -43.1),
        (8.25, -43.9),
    ]
    center = [(8.3, -43.0)]
    points = center + corners

    templates = get_templates_for_points(points, band)

    assert len(templates) == 7


def test_get_center_and_corners():
    image_path = os.path.join(os.path.dirname(__file__), "photometry_test_data", "RomanTDS", "images", "simple_model", "R062", "35083", "Roman_TDS_simple_model_R062_35083_8.fits.gz")
    points = get_center_and_corners(image_path)
    expected_df = pd.DataFrame(4, 24., 1, 23, 5, 16, 10, 293, 100, 234, 1, 14, 1, 1)
    pd.testing.assert_frame_equal(points, expected_df)
