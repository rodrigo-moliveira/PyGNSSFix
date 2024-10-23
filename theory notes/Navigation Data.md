GNSS pseudoranges are well known to be affected by instrumental biases. 
To cope with this problem, timing group delay (TGD) or differential 
code bias (DCB) parameters are commonly used for pseudorange corrections 
at the user ends. For the legacy GPS and GLONASS, it has been a common
practice to define clock offsets in both broadcast and precise
ephemeris products with respect to an ionosphere-free dual 
frequency combination of conventional reference signals
[L1/L2 P(Y)-code] (Montenbruck and Steigenberger 2013).
GPS/GLONASS applications using other signals or combined signals differing from the conventional reference signal
should apply TGD or DCB corrections, which are essential
for pseudorange-based positioning, timing and ionosphere
modeling (Wu et al. 2013; Montenbruck et al. 2014). Parameters designated as timing group delay (TGD) and inter-signal
corrections (ISC) in the navigation message are commonly
used to compensate the differential code biases in real time
for single-frequency users. Correction terms such as, TGD,
ISCL1C/A and ISCL2C, are initially defined to account for the
effect of space vehicle (SV) group delay differential between
L1 P(Y) and L2 P(Y), L1 P(Y) and L1 C/A, and between L1
P(Y) and L2 C, respectively (Tetewsky et al. 2009). In addition, more accurate differential code bias parameters designated as DCBs (such as DCBP1P2, DCBP1C1 and DCBP2C2)
are provided by GNSS communities to account for the same
delay as TGD and ISCs, particularly for the post-processing
applications (Schaer and Steigenberger 2006; Schaer 2008,
2012). The broadcast TGD values are referenced to an empirical absolute instrumental (satellite) bias, whereas DCB values are in a relative sense to reflect the differential hardware
(the satellite or receiver) delay between two different code
observations obtained on the same or two different frequencies (Li et al. 2012).



#### GPS Navigation Data ####
	The GPS TGD 
	The Timing Group Delay (TGD) is the bias difference between each GPS satellites P-code transmissions at the L1 and L2 frequencies. They are broadcast in the GPS navigation message so that single-frequency users can use them in conjunction with ionosphere delay estimates, such as the Klobuchar ionosphere model, to improve their position determination and to better derive UTC(USNO), which is the U.S. Naval Observatory’s realization of Coordinated Universal Time (UTC)
	By convention, the ionosphere-free linear combination (LC)
of P1 and P2 pseudorange is used for satellite clock esti-
mation, and the hardware delays are ignored during the
estimation process. As such, the GPS satellite clock error dts 
provided either in the broadcast navigation message or
the precise satellite clock product, contains a specific linear
combination of P1 and P2 satellite biases, specifically the
ionosphere-free LC

The L1 and L2 correction term, TGD, is initially calculated by the CS to account for the effect of SV group
delay differential between L1 P(Y) and L2 P(Y) based on measurements made by the SV contractor during
SV manufacture. The value of TGD for each SV may be subsequently updated to reflect the actual on-orbit
group delay differential. This correction term is only for the benefit of "single-frequency" (L1 C/A, L1 P(Y)
or L2 P(Y)) users; 

The value of TGD is not equal to the mean SV group delay differential, but is a measured value that
represents the mean group delay differential multiplied by 1/(1- gama). That is,
TGD = 1/(1-gama) * D(L1P-Y) - D(L2P-Y)
hardware biases for L1P(Y) and L2P(Y) signals respectively.
gama = (f_L1/f_L2)^2


P1 = rho + I + T + dt_r - (dt^s - d_P1) + err
P2 = rho + gama*I + T + dt_r - (dt^s - d_P2) + err
where:
	* rho is the true geometric range
	* err is the range error
	* I is ionosphere delay
	* gama is the frequency-dependent factor, gama=(f_1/f_2)^2
	* T is troposphere delay
	* dt^s is the satellite clock error
	* dt_r is the receiver clock error
	* d_P1, d_P2 are the instrumental (hardware) biases, which depend on the measured pseudorange (frequency)

By convention, the iionosphere-free linear combination (LC)
of P1 and P2 pseudorange is used for satellite clock estimation
and the hardware delays are ignored during the
estimation process. As such, the GPS satellite clock error
dt^eph provided either in the broadcast navigation message or
the precise satellite clock product, contains a specific linear
combination of P1 and P2 satellite biases, specifically the
ionosphere-free LC 

dt^eph = dt^s - (gama/(gama - 1) * D_P1 - 1/(gama - 1)*D_P2)

Therefore, an additional term to account for hardware
delay should be applied at the single-frequency user end
when using GPS satellite clock corrections. Otherwise, the
clock inconsistency will propagate into the positioning solution

P1 = rho + I + T + dt_r - dt^eph - 1/(gama - 1) * (D_P1 - D_P2) + err
P2 = rho + gama*I + T + dt_r - dt^eph - gama/(gama - 1) * (D_P1 - D_P2) + err

The dual-frequency iono-free combination is
P_IF(12) = rho + dt_r - dt^eph + err

 However, the instrumental
biases (DP1, DP2) cannot be determined in an absolute sense,
so the differential code biases are used instead
It is common to denote a specific difference of code
bias DCB(P1P2) as:
DCB(P1-P2) = DP1 − DP2

P1 = rho + I + T + dt_r - (dt^eph - TGD) + err
P2 = rho + gama*I + T + dt_r - (dt^eph - gama*TGD) + err


### Correction P1-C1

NOTE: the P-code is usually reported as 1W and 2W

For
maximum accuracy, coarse acquisition (C/A) code observations should first be converted to keep consistency with
P1/ P2 non-cross correlation types. A RINEX converter
utility, cc2noncc, is provided by IGS to easily make measurements consistent with P1/P2 data by applying satellitedependent P1–C1 bias corrections (Ray 2001).

from https://gssc.esa.int/navipedia/index.php/GPS_C1,_P1_and_P2_Codes_and_Receiver_Types
C1C + DCB(P1-C1) = C1P

C1P = rho + I + T + dt_r - (dt^s - d_P1) + err
C1C = rho + I + T + dt_r - (dt^s - d_C1) + err
DCB(P1-C1) = D_P1 - D_C1

-> C1P - C1C = D_P1 - D_C1 = DCB(P1-C1)
that is:
C1C + DCB(P1-C1) = C1P



### GALILEO
tenho de perceber se os BGDs / dt^eph são relativamente ao C1C-C5Q ou C1X-C5X, etc.
No gal não se põe este problema: não há DCBs entre 1C-1X ou 5Q-5X...