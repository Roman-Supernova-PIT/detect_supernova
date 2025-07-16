from detect_supernova.make_openuniverse_subtraction_pairs import get_image_info_for_ra_dec, get_templates_for_points


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
