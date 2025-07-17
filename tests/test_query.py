import os

import numpy as np
import pandas as pd

from detect_supernova.make_openuniverse_subtraction_pairs import (
    get_earliest_template_for_image,
    get_image_info_for_ra_dec,
    get_templates_for_points,
)
from detect_supernova.util import get_center_and_corners, make_data_records, read_data_records


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
    row = (4, 24.0, 1, 23, 5, 16, 10, 293, 100, 234)
    expected_df = pd.DataFrame.from_records([row], columns=expected_columns)

    points = get_center_and_corners(image_path)

    assert len(points) == 14

    # Make sure we have the columns we expect


#    assert len(set(expected_df.columns)) == len(set(expected_df.columns).intersection(set(points.columns)))


def test_get_earliest_template_for_image():
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

    points = get_center_and_corners(image_path)
    earliest_template = get_earliest_template_for_image(points)

    assert len(earliest_template) == 19


def test_get_image_info_for_ra_dec():
    expected_columns = (
        "boredec",
        "borera",
        "filter",
        "pointing",
        "sca",
        "ra",
        "dec",
        "ra_00",
        "dec_00",
        "ra_01",
        "dec_01",
        "ra_10",
        "dec_10",
        "ra_11",
        "dec_11",
    )
    expected_dtype = (
        np.float64,
        np.float64,
        str,
        np.int64,
        np.int64,
        np.float64,
        np.float64,
        np.float64,
        np.float64,
        np.float64,
        np.float64,
        np.float64,
        np.float64,
        np.float64,
        np.float64,
    )
    row = (1, 1, "R062", 1, 4, 24.0, 1, 23, 5, 16, 10, 293, 100, 234, 100)
    expected_df = pd.DataFrame.from_records([row], columns=expected_columns)

    ra, dec = 8.3, -42
    images = get_image_info_for_ra_dec(ra, dec)

    print(images.dtypes)
    assert len(images) == 484


def test_read_data_records_from_file():
    data_records_path = os.path.join(os.path.dirname(__file__), "test_ten_data_records.csv")
    data_records = read_data_records(data_records_path=data_records_path)

    assert len(data_records) == 10


def test_read_data_records_from_science_id_and_template_id():
    data_records = make_data_records(
        science_pointing=54670,
        science_sca=18,
        science_band="R062",
        template_pointing=26565,
        template_sca=18,
        template_band="R062",
        base_image_location=os.path.join(os.path.dirname(__file__), "photometry_test_data"),
    )

    assert data_records.template_pointing[0] == 26565
    assert len(data_records) == 1


def test_read_data_records_from_just_science_id():
    data_records = make_data_records(
        science_pointing=35083,
        science_sca=8,
        science_band="R062",
        base_image_location=os.path.join(os.path.dirname(__file__), "photometry_test_data"),
    )

    assert data_records.template_pointing[0] == 5044
    assert data_records.template_sca[0] == 8
    assert data_records.template_band[0] == "R062"
    assert len(data_records) == 1
