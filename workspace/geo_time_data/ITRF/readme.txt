Download the ITRF transformation files (with the Helmert Transformation coefficients) from:
https://epsg.io/

For instance, the transformation from ITRF93 to ITRF2000 is provided in:
https://epsg.io/9998

These rotations are implemented in this software to rotate the CSpice ephemerides (by default in ITRF93)
to the desired ITRF frame that is provided in the CODE IGS precise products (for example, in SP3 files).

The user may specify the desired ITRF frames in the configuration file:

```
    "cspice_frame": "ITRF93"
    "IGS_frame": "ITRF2020"
    "ITRF_rotation_file": "geo_time_data/ITRF/EPSG_9998.json"
```

Make sure that the ITRF rotation file is consistent with the provided ITRF frames.

NICE TO HAVE: this consistency could be checked automatically in the future by this software, and even download the
correct file from the EPSG website. This is future work.