import astropy.units as u

from coord_projection import two_direction_skymatch

MATCH_RADIUS = 0.4


def skymatch_and_join(
    left_table, right_table, left_skycoord, right_skycoord, match_radius=None
):
    match_radius = match_radius or MATCH_RADIUS

    left_table = left_table.copy().reset_index(drop=True)
    right_table = right_table.copy().reset_index(drop=True)

    matched_status, matched_id = two_direction_skymatch(
        left_skycoord, right_skycoord, radius=match_radius * u.arcsec
    )
    right_table = right_table.iloc[matched_id].copy().reset_index(drop=True)

    left_table["matched_status"] = matched_status
    joined_table = left_table.merge(right_table, left_index=True, right_index=True)

    return joined_table
