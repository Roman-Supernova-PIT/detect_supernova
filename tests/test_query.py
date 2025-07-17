import os

import numpy as np
import pandas as pd

from detect_supernova.make_openuniverse_subtraction_pairs import (
    get_image_info_for_ra_dec,
    get_templates_for_points,
    get_image_info_for_ra_dec,
)
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
    image_path = os.path.join(
        os.path.dirname(__file__),
        "photometry_test_data",
        "RomanTDS",
        "images",
        "simple_model",
        "R062",
        "35083",
        "Roman_TDS_simple_model_R062_35083_8.fits.gz",
    )
    expected_columns = ("ra", "dec", "ra_00", "dec_00", "ra_01", "dec_01", "ra_10", "dec_10", "ra_11", "dec_11")
    expected_dtype = (np.float64, np.float64, str, np.int64, np.int64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64)
    expected_df = pd.DataFrame.from_records(
        [(1, 1, "R062", 1, 4, 24.0, 1, 23, 5, 16, 10, 293, 100, 234)], columns=expected_columns
    )


    points = get_center_and_corners(image_path)

    assert len(points) == 1

    # Make sure we have the columns we expect
    assert len(set(expected_df.columns)) == len(set(expected_df.columns).intersection(set(points.columns)))


def test_get_image_info_for_ra_dec():
    expected_columns = ("boredec", "borera", "filter", "pointing", "sca", "ra", "dec", "ra_00", "dec_00", "ra_01", "dec_01", "ra_10", "dec_10", "ra_11", "dec_11")
    expected_dtype = (np.float64, np.float64, str, np.int64, np.int64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64)
    expected_df = pd.DataFrame.from_records(
        [(1, 1, "R062", 1, 4, 24.0, 1, 23, 5, 16, 10, 293, 100, 234)], columns=expected_columns
    )

    ra, dec = 8.3, -42
    images = get_image_info_for_ra_dec(ra, dec)

    print(images.dtype)
    assert len(images) == 484
