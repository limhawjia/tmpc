# This module converts Transverse Mercator projection coordinates (E, N) to their geographic equivalent (lat, long)
# Reference for the respective formulae is taken from:
# https://www.linz.govt.nz/data/geodetic-services/coordinate-conversion/projection-conversions/transverse-mercator-transformation-formulae

from math import pi as pi
from math import sin as sin
from math import radians as rad
from math import degrees as deg
from math import sqrt as sqrt
from math import tan as tan
from math import cos as cos


def sin2(x):
    return sin(x) * sin(x)


def sec(x):
    return 1 / cos(x)


# This class represents the origin of projection and is required by the Tmpc class to do any conversions.
class Origin:
    def __init__(self, latitude, longitude, false_easting, false_northing, central_meridian_scale_factor):
        self.latitude = latitude
        self.longitude = longitude
        self.false_easting = false_easting
        self.false_northing = false_northing
        self.central_meridian_scale_factor = central_meridian_scale_factor


# This class represents the reference spheroid and is required by the Tmpc class to do any conversions
class Spheroid:
    def __init__(self, semi_major_axis, ellipsoidal_flattening):
        self.semi_major_axis = semi_major_axis
        self.ellipsoidal_flattening = ellipsoidal_flattening


# This class performs the main conversion of the respective coordinates. Its constructor requires an Origin object
# representing the projection origin and a Spheroid object representing the reference spheroid for the particular
# Transverse Mercator projection coordinate system. Note that the values from the Origin and Spheroid objects
# are copied on instantiation
class Tmpc:
    def __init__(self, origin, spheroid):
        self.origin = origin
        self.spheroid = spheroid

    def convert_to_geographic(self, north, east):
        params = Parameters(self.origin, self.spheroid)

        # Calculate required terms
        north_prime = north - params.north_false
        m_prime = params.m0 + (north_prime / params.k)
        n = (params.a - params.b) / (params.a + params.b)
        n2 = n * n
        n3 = n2 * n
        n4 = n2 * n2
        g = params.a * (1 - n) * (1 - n2) * (1 + 9 * n2 / 4 + 255 * n4 / 64) * (pi / 180)
        sigma = (m_prime * pi) / (180 * g)
        lat_prime_term1 = (3 * n / 2 - 27 * n3 / 32) * sin(2 * sigma)
        lat_prime_term2 = (21 * n2 / 16 - 55 * n4 / 32) * sin(4 * sigma)
        lat_prime_term3 = 151 * n3 / 96 * sin(6 * sigma)
        lat_prime_term4 = 1097 * n4 / 512 * sin(8 * sigma)
        lat_prime = sigma + lat_prime_term1 + lat_prime_term2 + lat_prime_term3 + lat_prime_term4

        rho_prime = params.a * (1 - params.e2) / pow(1 - params.e2 * sin2(lat_prime), 3 / 2)
        v_prime = params.a / sqrt(1 - params.e2 * sin2(lat_prime))
        psi_prime = v_prime / rho_prime
        psi_prime2 = psi_prime * psi_prime
        psi_prime3 = psi_prime * psi_prime * psi_prime
        psi_prime4 = psi_prime2 * psi_prime2
        t_prime = tan(lat_prime)
        t_prime2 = t_prime * t_prime
        t_prime4 = t_prime2 * t_prime2
        t_prime6 = t_prime2 * t_prime2 * t_prime2
        east_prime = east - params.east_false
        x = east_prime / (params.k * v_prime)
        x3 = x * x * x
        x5 = x3 * x * x
        x7 = x5 * x * x

        # Get latitude
        lat_factor = t_prime / (params.k * rho_prime)
        lat_term1 = lat_factor * (east_prime * x) / 2
        lat_term2 = lat_factor * (east_prime * x3) / 24 * (
                -4 * psi_prime + 9 * psi_prime * (1 - t_prime2) + 12 * t_prime2)
        lat_term3 = lat_factor * (east_prime * x5) / 720 * \
                    (8 * psi_prime4 * (11 - 24 * t_prime4) -
                     12 * psi_prime3 * (21 - 71 * t_prime2) +
                     15 * psi_prime2 * (15 - 98 * t_prime2 + 15 * t_prime4) +
                     180 * psi_prime * (5 * t_prime2 - 3 * t_prime4) +
                     360 * t_prime4)
        lat_term4 = lat_factor * (east_prime * x7) / 40320 * (
                1385 + 3633 * t_prime2 + 4095 * t_prime4 + 1575 * t_prime6)
        lat_radians = lat_prime - lat_term1 + lat_term2 - lat_term3 + lat_term4
        lat = deg(lat_radians)

        # Get longitude
        long_term1 = x * sec(lat_prime)
        long_term2 = x3 * sec(lat_prime) / 6 * (psi_prime + 2 * t_prime2)
        long_term3 = x5 * sec(lat_prime) / 120 * \
                     (-4 * psi_prime3 * (1 - 6 * t_prime2) +
                      psi_prime2 * (9 - 68 * t_prime2) +
                      72 * psi_prime * t_prime2 +
                      24 * t_prime4)
        long_term4 = x7 * sec(lat_prime) / 5040 * \
                     (61 + 622 * t_prime2 + 1320 * t_prime4 + 720 * t_prime6)
        long_radians = rad(params.long_origin) + long_term1 - long_term2 + long_term3 - long_term4
        long = deg(long_radians)

        return round(lat, 10), round(long, 10)

    def convert_to_projection(self):
        raise NotImplementedError


# This class represents the entire set of common parameters needed conversion between the two coordinate systems. They
# include values copied form the origin of projection and reference spheroid upon instantiation along with other
# derived parameters.
class Parameters:
    def __init__(self, origin, spheroid):
        self.a = spheroid.semi_major_axis
        self.f = spheroid.ellipsoidal_flattening
        self.lat_origin = origin.latitude
        self.long_origin = origin.longitude
        self.north_false = origin.false_northing
        self.east_false = origin.false_easting
        self.k = origin.central_meridian_scale_factor

        self.b = self.a * (1 - self.f)
        self.e2 = 2 * self.f - self.f * self.f
        self.e4 = self.e2 * self.e2
        self.e6 = self.e2 * self.e2 * self.e2

        A0 = 1 - self.e2 / 4 - 3 * self.e4 / 64 - 5 * self.e6 / 256
        A2 = 3 / 8 * (self.e2 + self.e4 / 4 + 15 * self.e6 / 128)
        A4 = 15 / 256 * (self.e4 + 3 * self.e6 / 4)
        A6 = 35 * self.e6 / 3072

        self.m0 = self.a * (A0 * rad(self.lat_origin) -
                            A2 * sin(2 * rad(self.lat_origin)) +
                            A4 * sin(4 * rad(self.lat_origin)) -
                            A6 * sin(6 * rad(self.lat_origin)))
