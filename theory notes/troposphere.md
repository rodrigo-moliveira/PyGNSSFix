Tropospheric Delay

The Zenith Total Delay (ZTD) is:
    ZTD = ZHD + ZWD
where
    * ZHD -> Zenith Hydrostatic Delay
    * ZWD -> Zenith Wet Delay

With an appropriate mapping function and using the satellite elevation angle theta as input, the
ZTD can be resolved into the Slant Tropospheric Delay (STD) as:
    STD = ZHD * MF_h(theta) + ZWD * MF_w
where MF_h and MF_w represent hydrostatic and wet mapping functions respectively



# USAGE OF GPT TROPO MODELS

1) Use GPT2w + VMF1
Get file gpt2_1wA.grd and read grid with read_grid script
Use gpt2_1w.m to get p, T, dT, ah, aw...
Use vmf1_ht.m to get vmf1h and vmf1w mapping coefficients
Apply (2) to get zhd and (3) to get zwd
Apply (1) to get ztd = zhd*vmf1h + zwd*vmf1w

De forma semelhante, usar o gpt3_5.grd, aplicar o modelo com gpt3_5.fast.py
Use vmf3_ht.py to get vmf3h and vmf3w (numa primeira fase posso usar o VMF1, que já está me python)

Também posso usar a GMF como mapping function



# Compute ZWD and ZHD
ZHD = 0.0022768 * p / (1 - 0.00266 * cos(2*lat) - 0.28e-6 * height)  #Eq. (2)
where:
    * p is pressure in mbar (1hPa = 1mbar)
    * lat is latitude in rad
    * height is in m

    ZTD = 1e-6 * (k2_ + k3/Tm) * (Rd * e) / (gm * (lambda + 1))
where:
    * k2 constant equal to 16.52K/mbar
    * k3 constant equal to 3776E5 K^2/mbar
    * Tm mean temperature weighted with the water vapor in [K]
    * Rd specific gas constant of dry air = 287.0464 JK^-1kg^-1
    * e water vapor pressure in mbar
    * gm local gravity = 9.80665 ms^-2
    * lambda water decrease factor [dimensionless]