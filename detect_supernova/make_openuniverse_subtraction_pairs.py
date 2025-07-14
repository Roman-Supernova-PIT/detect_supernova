"""Find science, templates pairs for example RA, Dec for OpenUniverse2024 images"""

import argparse
import requests

import pandas as pd


def get_image_info_for_ra_dec(ra, dec, band=None):
    server_url = "https://roman-desc-simdex.lbl.gov"
    req = requests.Session()
    json = {"containing": [ra, dec]}
    if band is not None:
        json["filter"] = band
    result = req.post(f"{server_url}/findromanimages", json=json)
    if result.status_code != 200:
        raise RuntimeError(f"Got status code {result.status_code}\n{result.text}")
    df = pd.DataFrame(result.json())
    return df


def get_templates_for_points(points, band, min_points=3):
    """Returns all images in the same bandpass that overlap at least min_points
    out of the list of points passed in.

    Parameters
    ----------
    points : List of tuples of (ra, dec) points
    band : str
    min_points : int, number of required points images must cover

    Returns
    -------
    images : List of DataFrames of image info from Roman-DESC-simdex

    Notes
    -----
    This uses the roman-desc-simdex NERSC Spin server.
    The DB server query for image info as a function of RA, Dec sometimes fails.
    If it does then it will raise a RuntimeError, this function will let that pass through.
    """
    matches = []
    for i, (ra, dec) in enumerate(points):
        matching_images = get_image_info_for_ra_dec(ra, dec, band=band)
        matches.append(matching_images)

    matches = pd.concat(matches)
    # From
    #  https://stackoverflow.com/questions/35584085/how-to-count-duplicate-rows-in-pandas-dataframe
    matches = matches.groupby(matches.columns.tolist()).size().reset_index().rename(columns={0: "counts"})
    good_matches = matches.loc[matches.counts >= min_points]

    return good_matches


def get_templates_for_image(im, min_points=3):
    """Return a list of matching images that could be used as templates.

    Returns all images in the same bandpass that overlap at least min_points
    out of the 5 points of the center + corners of the images

    Parameters
    ----------
    images: Object with data attributes
    ("ra", "dec", "ra_00", "dec_00", "ra_01", "dec_01", "ra_10", "dec_10", "ra_11", "dec_11")
    and get method for "filter"
    min_points: int

    Returns
    -------
    images : list of (pointing, sca, band) tuples of overlapping images
    """
    corners = [
        (im.ra_00, im.dec_00),
        (im.ra_01, im.dec_01),
        (im.ra_10, im.dec_10),
        (im.ra_11, im.dec_11),
    ]
    center = [(im.ra, im.dec)]
    points = center + corners

    # band.filter would be a method so we can't use data attribute of same name
    # and instead use `get` to access
    band = im.get("filter")

    return get_templates_for_points(points, band=band, min_points=min_points)


def get_earliest_template_for_image(image, **kwargs):
    """
    Parameters
    ----------
    image : DataFrame of image info from Roman-DESC-simdex
    """
    templates = get_templates_for_image(image, **kwargs)
    # Get earliest MJD
    earliest_template = templates.iloc[templates.mjd.argsort()].iloc[0]

    return earliest_template


def run(
    band="R062",
    transient_type="astrophysical",
    ra=[8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0],
    dec=[-43.2, -43.2, -43.2, -43.2, -43.2, -43.1, -43.1, -43.1, -43.1, -43.1],
):
    """Get templates for 10 different pointings

    I don't know how to query the roman-desc-simdex for a given pointing, sca

    So I'm just going to do the silly thing and ask for
    10 different RA, DEC positions and use that image info.

    From
    https://arxiv.org/pdf/2501.05632
    The RomanWAS MJD range is 61444 -- 63270
    The RomanTDS MJD range is 62000 -- 63563

    "astrophysical" are simulated based on 10 different types of astrophysical transients
    from 61444, 63269.
    The peak MJD range was from 61374 to 63299 so that the simulated range is complete

    "ranmag" is a flat AB spectrum -22 < mag < -17 from 0 < z < 2.
    meant for easy testing of subtractions
    The ranmag are just the last 14 days of the simulated range.
    """
    mjd_range = {"astrophysical": (61444, 63269), "ranmag": (63549, 63563)}
    mjd_min, mjd_max = mjd_range[transient_type]

    science_images = []
    for r, d in zip(ra, dec):
        image_info = get_image_info_for_ra_dec(r, d, band=band)
        # Take the last one by MJD that's within the MJD transient range
        # of the OpenUniverse sims
        image_info = image_info.loc[(mjd_min < image_info.mjd) & (image_info.mjd < mjd_max)]
        science_image = image_info.iloc[image_info.mjd.argsort()].iloc[-1]
        science_images.append(science_image)

    print(science_images)
    # Now we have 10 images with the information we need to get center and corners.
    # Now we're going to pick the earliest template for each science image
    #  that overlaps at least 3 (corners, center).
    template_images = []
    for im in science_images:
        template = get_earliest_template_for_image(im)
        template_images.append(template)

    head = "science_band,science_pointing,science_sca,template_band,template_pointing,template_sca"
    print(head)
    for s, t in zip(science_images, template_images):
        print(s.get("filter"), s.pointing, s.sca, t.get("filter"), t.pointing, t.sca)

    return science_images, template_images


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--band", "-b", type=str, default="R062")
    parser.add_argument(
        "--transient_type",
        "-t",
        type=str,
        choices=["astrophysical", "ranmag"],
        default="astrophysical",
        help="Choose simulated transient class to specify MJD range based on different type of injected transient.",
    )
    args = parser.parse_args()

    run(**args.__dict__)
