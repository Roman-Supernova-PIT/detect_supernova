"""
Find templates for given images
"""

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
    """
    Returns all images in the same bandpass that overlap at least min_points
    out of the list of points passed in.

    Parameters
    ----------
    points : List of tuples of (ra, dec) points
    band : str
    min_points : int, number of required points images must cover

    Returns
    -------
    images : list of (pointing, sca, band) tuples of overlapping images
    """
    matches = []
    for i, (ra, dec) in enumerate(points):
        matching_images = get_image_info_for_ra_dec(ra, dec, band=band)
        matches.append(matching_images)

    matches = pd.concat(matches)
    # From
    #  https://stackoverflow.com/questions/35584085/how-to-count-duplicate-rows-in-pandas-dataframe
    matches = (
        matches.groupby(matches.columns.tolist())
        .size()
        .reset_index()
        .rename(columns={0: "counts"})
    )
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


def run():
    """Get templates for 10 different pointings

    I don't know how to easily query the roman-desc-simdex for a given pointing, sca

    So I'm just going to do the silly thing and ask for 10 different RA, DEC positions
    and use that image info.
    """
    ra = [8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0]
    dec = [-43.2, -43.2, -43.2, -43.2, -43.2, -43.1, -43.1, -43.1, -43.1, -43.1]
    band = "R062"

    science_images = []
    for r, d in zip(ra, dec):
        image_info = get_image_info_for_ra_dec(r, d, band=band)
        # Take the last one by MJD
        science_image = image_info.iloc[image_info.mjd.argsort()].iloc[0]
        science_images.append(science_image)

    # Now we have 10 images with the information we need to get center and corners.
    # Now we're going to We're just going to pick one for each and have our list.
    template_images = []
    for im in science_images:
        possible_templates = get_templates_for_image(im)
        # Get earliest MJD
        template = possible_templates.iloc[possible_templates.mjd.argsort()].iloc[0]
        template_images.append(template)

    head = "science_band,science_pointing,science_sca,template_band,template_pointing,template_sca"
    print(head)
    for s, t in zip(science_images, template_images):
        print(s.get("filter"), s.pointing, s.sca, t.get("filter"), t.pointing, t.sca)

    return science_images, template_images


if __name__ == "__main__":
    run()
