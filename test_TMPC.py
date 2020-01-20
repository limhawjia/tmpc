from tmpc.tmpc import *

svy21_origin = Origin(1.3666666666, 103.83333333, 28001.6420, 38744.5720, 1.0000)
svy21_spheroid = Spheroid(6378137.0000, 1 / 298.257223563)

converter = Tmpc(svy21_origin, svy21_spheroid)

# SVY21 coordinates for Arab Street
easting = 30770.0173
northing = 31622.8216

latitude = 1.302253804056663
longitude = 103.85820845283774

lat, long = converter.convert_to_geographic(northing, easting)
print(lat)
print(long)
