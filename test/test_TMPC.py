from tmpc import *


def within_threshold(threshold, expected, actual):
    return expected - threshold < actual < expected + threshold


def test_convert_to_geographic():
    geographic_threshold = 0.00001
    projection_threshold = 1

    # Singapore's svy21 system's origin of projection and reference spheroid
    svy21_origin = Origin(1.3666666666, 103.83333333, 28001.6420, 38744.5720, 1.0000)
    svy21_spheroid = Spheroid(6378137.0000, 1 / 298.257223563)

    converter = Tmpc(svy21_origin, svy21_spheroid)

    # Arab Street, Singapore
    east1 = 30770.0173
    north1 = 31622.8216
    lat1 = 1.302253804056663
    long1 = 103.85820845283774

    result_lat1, result_long1 = converter.convert_to_geographic(north1, east1)

    assert within_threshold(geographic_threshold, lat1, result_lat1)
    assert within_threshold(geographic_threshold, long1, result_long1)

    result_north1, result_east1 = converter.convert_to_projection(lat1, long1)

    assert within_threshold(projection_threshold, east1, result_east1)
    assert within_threshold(projection_threshold, north1, result_north1)
