# TMPC

This module provides the class `Tmpc` which converts coordinates between the geographic coordinates system (lat, long) and the Transverse Mercator Projection coordinate system (northing, easting).

The Transverse Mercator Projection coordinate system is a coordinate system that is widely used in national and international mapping systems. When paired with a suitable geodetic statum -- an origin of projection and a reference spheroid, it delivers high accuracies in zones less than a few degrees in the East-West extent.

# Credits

Basic information : https://en.wikipedia.org/wiki/Transverse_Mercator_projection

In-depth information : https://gisgeography.com/utm-universal-transverse-mercator-projection/

Transfomration formulae : https://www.linz.govt.nz/data/geodetic-services/coordinate-conversion/projection-conversions/transverse-mercator-transformation-formulae

# Usage

```python
from tmpc import *

origin = Origin(1.3666666666, 103.83333333, 28001.6420, 38744.5720, 1.0000)
spheroid = Spheroid(6378137.0000, 1 / 298.257223563)
converter = Tmpc(origin, spheroid)

northing = 31622.8216
easting = 30770.0173

lat, long = converter.convert_to_geographic(northing, easting)
```

`Tmpc` is the main class used to convert coordinates from the two system. An origin of projection, `Origin` and reference spheroid, `Spheroid` is required for instantiation.

The constructors for `Origin` and `Spheroid` are as follows:

Origin:

```python
def __init__(self, latitude, longitude, false_easting, false_northing, central_meridian_scale_factor)
```

Spheroid:

```python
def __init__(self, semi_major_axis, ellipsoidal_flattening)
```

To find out more about the different parameters, visit the links above to learn about the Transverse Mercator Projection coordinate system. Be sure to choose values that corresponds the the specific system that you are dealing with. 

For example, the data for Singapore's SVY21 coordinate system can be found here:

https://www.sla.gov.sg/sirent/CoordinateSystems.aspx

# Disclaimer
The conversion results in some loss of precision. This module has also not been tested extensively enough to gurantee 100% accuracy.
