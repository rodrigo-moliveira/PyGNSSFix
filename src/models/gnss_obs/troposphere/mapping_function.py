import numpy as np


class VMF1:
    @staticmethod
    def compute(ah=None, aw=None, dmjd=None, dlat=None, ht=None, zd=None):
        # =============================================================================
        #     !!! This is the version with height correction !!!
        #     !!! It has to be used with the grid !!!
        #
        #     This subroutine determines the VMF1 (Vienna Mapping Functions 1)
        #     Reference: Boehm, J., B. Werl, H. Schuh (2006),
        #     Troposphere mapping functions for GPS and very long baseline interferometry
        #     from European Centre for Medium-Range Weather Forecasts operational analysis data,
        #     J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.
        #
        #     Please mind that the coefficients in this paper are wrong. The corrected version of
        #     the paper can be found at:
        #     http://vmf.geo.tuwien.ac.at/documentation/Boehm%20et%20al.,%202006a%20%20%20(VMF1).pdf
        #
        #     input data
        #     ----------
        #     ah:   hydrostatic coefficient a (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
        #     aw:   wet coefficient a         (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
        #     dmjd: modified julian date
        #     dlat: ellipsoidal latitude in radians
        #     ht:   ellipsoidal height in meter
        #     zd:   zenith distance in radians
        #
        #     output data
        #     -----------
        #     vmf1h: hydrostatic mapping function
        #     vmf1w: wet mapping function
        #
        #     Johannes Boehm, 2005 October 2
        #      Rev 2011 July 21: latitude -> ellipsoidal latitude

        #     implicit double precision (a-h,o-z)

        #     pi = 3.14159265359d0

        #     reference day is 28 January
        #     this is taken from Niell (1996) to be consistent
        #     File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
        # =============================================================================
        doy = dmjd - 44239.0 + 1 - 28
        bh = 0.0029
        c0h = 0.062
        if dlat < 0:
            phh = np.pi
            c11h = 0.007
            c10h = 0.002
        else:
            phh = 0.0
            c11h = 0.005
            c10h = 0.001

        ch = c0h + ((np.cos(doy / 365.25 * 2.0 * np.pi + phh) + 1.0) * c11h / 2.0 + c10h) * (1.0 - np.cos(dlat))
        sine = np.sin(np.pi / 2.0 - zd)
        beta = bh / (sine + ch)
        gamma = ah / (sine + beta)
        topcon = (1.0 + ah / (1.0 + bh / (1.0 + ch)))
        vmf1h = topcon / (sine + gamma)
        # C  height correction for hydrotatic part [Niell, 1996]
        a_ht = 2.53e-05
        b_ht = 0.00549
        c_ht = 0.00114
        hs_km = ht / 1000.0
        beta = b_ht / (sine + c_ht)
        gamma = a_ht / (sine + beta)
        topcon = (1.0 + a_ht / (1.0 + b_ht / (1.0 + c_ht)))
        ht_corr_coef = 1.0 / sine - topcon / (sine + gamma)
        ht_corr = np.multiply(ht_corr_coef, hs_km)
        vmf1h = vmf1h + ht_corr
        bw = 0.00146
        cw = 0.04391
        beta = bw / (sine + cw)
        gamma = aw / (sine + beta)
        topcon = (1.0 + aw / (1.0 + bw / (1.0 + cw)))
        vmf1w = topcon / (sine + gamma)
        return vmf1h, vmf1w


class VMF3:
    @staticmethod
    def vmf3_ht(ah, aw, mjd, lat, lon, h_ell, zd):
        """
            !!! This is the version with height correction !!!
            !!! It has to be used with the grid !!!
        (c) Department of Geodesy and Geoinformation, Vienna University of
        Technology, 2016
        This subroutine determines the VMF3 hydrostatic and wet mapping factors.
        The a coefficients have to be inserted from discrete data, while the b
        and c coefficients are of empirical nature containing a geographical
        and temporal dependence, represented in spherical harmonics. The
        spherical harmonics coefficients are developed to degree and order 12 and
        are based on a 5�x5� grid containing ray-tracing data from 2001-2010.
        All input quantities have to be scalars!
        INPUT:
            ah: hydrostatic mf coefficient a (http://vmf.geo.tuwien.ac.at/trop
            aw: wet mf coefficient a (http://vmf.geo.tuwien.ac.at/trop_product
            mjd: modified Julian date
            lat: latitude (radians)
            lon: longitude (radians)
            h_ell: ellipsoidal height (m)
            zd: zenith distance (radians)
        OUTPUT:
            mfh: hydrostatic mapping factor
            mfw: wet mapping factor
        """
        # convert mjd to doy

        hour = int(np.floor((mjd - int(np.floor(mjd))) * 24))
        minu = int(np.floor((((mjd - int(np.floor(mjd))) * 24) - hour) * 60))
        sec = (((((mjd - int(np.floor(mjd))) * 24) - hour) * 60) - minu) * 60

        # change secs, min hour whose sec==60

        if sec == 60:
            minu = minu + 1
            sec = 0

        if minu == 60:
            hour = hour + 1
            minu = 0
        # calc jd (yet wrong for hour==24)
        jd = mjd + 2400000.5
        # if hr==24, correct jd and set hour==0
        if hour == 24:
            jd = jd + 1
            hour = 0

        # integer julian date
        jd_int = int(np.floor(jd + 0.5))
        aa = jd_int + 32044
        bb = int(np.floor((4 * aa + 3) / 146097))
        cc = aa - int(np.floor((bb * 146097) / 4))
        dd = int(np.floor((4 * cc + 3) / 1461))
        ee = cc - int(np.floor((1461 * dd) / 4))
        mm = int(np.floor((5 * ee + 2) / 153))
        day = ee - int(np.floor((153 * mm + 2) / 5)) + 1
        month = mm + 3 - 12 * int(np.floor(mm / 10))
        year = bb * 100 + dd - 4800 + int(np.floor(mm / 10))

        # first check if the specified year is leap year or not (logical output)
        leapYear = (np.logical_or((np.mod(year, 4) == np.logical_and(0, np.mod(year, 100)) != 0),
                                  np.mod(year, 400)) == 0)
        days = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        doy = sum(days[range(0, month - 1)]) + day
        if leapYear == 1 and month > 2:
            doy = doy + 1

        doy = doy + mjd - int(np.floor(mjd))

        # https://chat.openai.com/c/57e0c54e-0fb2-4672-9882-aeef45186666
