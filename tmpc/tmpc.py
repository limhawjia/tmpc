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
    return 1.0 / cos(x)


def cos2(x):
    return cos(x) * cos(x)


def cos3(x):
    return cos(x) * cos(x) * cos(x)


def cos4(x):
    return cos2(x) * cos2(x)


def cos5(x):
    return cos3(x) * cos(x) * cos(x)


def cos6(x):
    return cos3(x) * cos3(x)


def cos7(x):
    return cos5(x) * cos(x) * cos(x)


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
        sigma = (m_prime * pi) / (180.0 * g)
        lat_prime_term1 = (3.0 * n / 2.0 - 27.0 * n3 / 32.0) * sin(2.0 * sigma)
        lat_prime_term2 = (21.0 * n2 / 16.0 - 55.0 * n4 / 32.0) * sin(4.0 * sigma)
        lat_prime_term3 = 151.0 * n3 / 96.0 * sin(6.0 * sigma)
        lat_prime_term4 = 1097.0 * n4 / 512.0 * sin(8.0 * sigma)
        lat_prime = sigma + lat_prime_term1 + lat_prime_term2 + lat_prime_term3 + lat_prime_term4

        rho_prime = params.a * (1.0 - params.e2) / pow(1.0 - params.e2 * sin2(lat_prime), 3.0 / 2.0)
        v_prime = params.a / sqrt(1.0 - params.e2 * sin2(lat_prime))
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
        lat_term1 = lat_factor * (east_prime * x) / 2.0
        lat_term2 = lat_factor * (east_prime * x3) / 24.0 * (
                -4.0 * psi_prime + 9.0 * psi_prime * (1.0 - t_prime2) + 12.0 * t_prime2)
        lat_term3 = lat_factor * (east_prime * x5) / 720.0 * \
                    (8.0 * psi_prime4 * (11.0 - 24.0 * t_prime4) -
                     12.0 * psi_prime3 * (21.0 - 71.0 * t_prime2) +
                     15.0 * psi_prime2 * (15.0 - 98.0 * t_prime2 + 15.0 * t_prime4) +
                     180.0 * psi_prime * (5.0 * t_prime2 - 3.0 * t_prime4) +
                     360.0 * t_prime4)
        lat_term4 = lat_factor * (east_prime * x7) / 40320.0 * (
                1385.0 + 3633.0 * t_prime2 + 4095.0 * t_prime4 + 1575.0 * t_prime6)
        lat_radians = lat_prime - lat_term1 + lat_term2 - lat_term3 + lat_term4
        lat = deg(lat_radians)

        # Get longitude
        long_term1 = x * sec(lat_prime)
        long_term2 = x3 * sec(lat_prime) / 6.0 * (psi_prime + 2.0 * t_prime2)
        long_term3 = x5 * sec(lat_prime) / 120.0 * \
                     (-4.0 * psi_prime3 * (1.0 - 6.0 * t_prime2) +
                      psi_prime2 * (9.0 - 68.0 * t_prime2) +
                      72.0 * psi_prime * t_prime2 +
                      24.0 * t_prime4)
        long_term4 = x7 * sec(lat_prime) / 5040.0 * \
                     (61.0 + 622.0 * t_prime2 + 1320.0 * t_prime4 + 720.0 * t_prime6)
        long_radians = rad(params.long_origin) + long_term1 - long_term2 + long_term3 - long_term4
        long = deg(long_radians)

        return round(lat, 10), round(long, 10)

    def convert_to_projection(self, lat, long):
        params = Parameters(self.origin, self.spheroid)

        m = params.a * (params.A0 * rad(lat) -
                        params.A2 * sin(2 * rad(lat)) +
                        params.A4 * sin(4 * rad(lat)) -
                        params.A6 * sin(6 * rad(lat)))
        rho = params.a * (1 - params.e2) / pow(1 - params.e2 * sin2(rad(lat)), 3.0 / 2.0)
        v = params.a / sqrt(1 - params.e2 * sin2(rad(lat)))
        psi = v / rho
        psi2 = psi * psi
        psi3 = psi * psi * psi
        psi4 = psi2 * psi2
        t = tan(rad(lat))
        t2 = t * t
        t4 = t2 * t2
        t6 = t4 * t2
        omega = rad(long - params.long_origin)
        omega2 = omega * omega
        omega4 = omega2 * omega2
        omega6 = omega2 * omega4
        omega8 = omega4 * omega4

        north_term1 = omega2 / 2.0 * v * sin(rad(lat)) * cos(rad(lat))
        north_term2 = omega4 / 24.0 * v * sin(rad(lat)) * cos3(rad(lat)) * (4 * psi2 + psi - t2)
        north_term3 = omega6 / 720 * v * sin(rad(lat)) * cos5(rad(lat)) * \
                      (8 * psi4 * (11 - 24 * t2) -
                       28 * psi3 * (1 - 6 * t2) +
                       psi2 * (1 - 32 * t2) -
                       psi * 2 * t2 +
                       t4)
        north_term4 = omega8 / 40320 * v * sin(rad(lat)) * cos7(rad(lat)) * \
                      (1385 - 3111 * t2 + 543 * t4 - t6)
        north = params.north_false + params.k * (m - params.m0 + north_term1 + north_term2 + north_term3 + north_term4)

        east_term1 = omega2 / 6 * cos2(rad(lat)) * (psi - t2)
        east_term2 = omega4 / 120 * cos4(rad(lat)) * \
                     (4 * psi3 * (1 - 6 * t2) +
                      psi2 * (1 + 8 * t2) -
                      psi * 2 * t2 +
                      t4)
        east_term3 = omega6 / 5040 * cos6(rad(lat)) * (61 - 479 * t2 + 179 * t4 - t6)
        east = params.east_false + params.k * v * omega * cos(rad(lat)) * (1 + east_term1 + east_term2 + east_term3)

        return round(north, 10), round(east, 10)


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

        self.b = self.a * (1.0 - self.f)
        self.e2 = 2.0 * self.f - self.f * self.f
        self.e4 = self.e2 * self.e2
        self.e6 = self.e2 * self.e2 * self.e2

        self.A0 = 1.0 - self.e2 / 4.0 - 3.0 * self.e4 / 64.0 - 5.0 * self.e6 / 256.0
        self.A2 = 3.0 / 8.0 * (self.e2 + self.e4 / 4.0 + 15.0 * self.e6 / 128.0)
        self.A4 = 15.0 / 256.0 * (self.e4 + 3.0 * self.e6 / 4.0)
        self.A6 = 35.0 * self.e6 / 3072.0

        self.m0 = self.a * (self.A0 * rad(self.lat_origin) -
                            self.A2 * sin(2 * rad(self.lat_origin)) +
                            self.A4 * sin(4 * rad(self.lat_origin)) -
                            self.A6 * sin(6 * rad(self.lat_origin)))
