from tmpc import *
import sys

threshold = 0.00001


def within_threshold(expected, actual):
    return expected - threshold < actual < expected + threshold


def test_convert_to_geographic():
    #
    svy21_origin = Origin(1.3666666666, 103.83333333, 28001.6420, 38744.5720, 1.0000)
    svy21_spheroid = Spheroid(6378137.0000, 1 / 298.257223563)

    converter = Tmpc(svy21_origin, svy21_spheroid)

    # Arab Street, Singapore
    east1 = 30770.0173
    north1 = 31622.8216
    lat1 = 1.302253804056663
    long1 = 103.85820845283774

    actual_lat1, actual_long1 = converter.convert_to_geographic(north1, east1)

    assert within_threshold(lat1, actual_lat1)
    assert within_threshold(long1, actual_long1)



