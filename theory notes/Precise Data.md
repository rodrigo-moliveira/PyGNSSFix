# IGS Clock Products (pag 1007 of handbook):
The GPS clock parameters provided by the IGS
refer (by convention) to the ionosphere-free linear combination of the L1 P(Y) and L2 P(Y) signals. However
several geodetic GPS receivers only track the C/A code
on L1. Therefore, DCBs between the different types of
signals have to be considered (Sect. 19.6.1). One can
either estimate the DCB parameters, as done by the
CODE AC, or one can correct for the DCBs, for example, with the cc2noncc tool utilizing the CODE DCB
estimates [34.79]. Users of IGS clock products must ensure consistency between the clock parameters on the
one hand and the observation types on the other hand
by applying the corresponding DCBs. Further sources
for DCBs are the Time Group Delay (TGD) parameters of the GPS navigation message (Sect. 7.4.3), the
Inter Signal Corrections (ISCs) of the GPS Civil Navigation Message CNAV [34.80], and the IGS MGEX
DCB product [34.81].

For GAL it is E1/E5a

### DCB Files


Differential Code Biases (DCBs) are the systematic errors, or biases, between two GNSS code observations at the same or different frequencies. DCBs are required for code-based positioning of GNSS receivers, extracting ionosphere total electron content (TEC), and other applications. Proper knowledge of DCBs is crucial to many navigation applications but also non-navigation applications such as ionospheric analysis and time transfer. With all of the new signals offered by modernized and new GNSS constellations, analysts now require a comprehensive multi-GNSS DCB product.

As part of the IGS Multi-GNSS Experiment (MGEX) two different DCB products are determined by two analysis groups: Institute of Geodesy and Geophysics (IGG) of the Chinese Academy of Sciences (CAS) in Wuhan and the Deutsche Forschungsanstalt für Luftund Raumfahrt (DLR) in Germany. These files are archived at the CDDIS:

https://cddis.nasa.gov/archive/gnss/products/bias/

The satellite and station biases computed by the Institute of Geodesy and Geophysics (IGG) of the Chinese Academy of Sciences (CAS) in Wuhan are generated on a daily basis with a latency of 2-3 days and made available in daily BSX files. The constellations and signals analyzed by the IGG CAS are



The DLR product is updated on a 3-monthly basis and comprises weekly averages of the satellite biases. In addition a full set of daily satellite and station biases is provided for reference purposes and trend analyses. The constellations and signals analyzed by the DLR are:
source https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/gnss_differential_code_bias_product.html


 GPS C1C observations have to be corrected according to                         
   Obs_C1C - Bias_C1C + Bias_C1W                                                
 to make them consistent to C1W. This corresponds to the cc2noncc correction:   
   Obs_C1C + DCB                                                                
 with the classical GPS P1-C1 DCB = Bias_C1W - Bias_C1C.                        
 The correction is similar for old-fashioned GPS cross-correlation (cc)         
 observations C2D = C1C + (C2W - C1W) to make them consistent to C2W:           
   Obs_C2D + DCB = Obs_C2D + (Bias-C1W - Bias_C1C).  







Within the Multi-GNSS Experiment (MGEX) launched by the International GNSS Service (IGS) (Dow et al. 2009), over 50 stations are now
equipped with new generation multi-system receivers that
can record BDS signals worldwide (Rizos et al. 2013; Dach
et al. 2014). These offer a basis for the determination
of BDS DCBs. Currently both BDS satellite and receiver
biases from weekly averages of daily DCBs are provided
in the annual files MGEXyyyy.bsx and MGEXyyyy_all.bsx,
where yyyy indicates the four-digit year. Figure 1 shows
parts of the format of DCBs extracted from the annual file
“MGEX2014.bsx” (available at: ftp://cddis.gsfc.nasa.gov/
pub/gps/products/mgex/dcb/)




BIAS | Bias type identifier. | 1X,A4 |
| | Available types are: | |
| | ’DSB ’: Differential Signal | |
| | Bias (DSB); | |
| | ’ISB ’: Ionosphere-free (linear | |
| | combination) Signal | |
| | Bias (ISB); | |
| | ’OSB ’: Observable-specific | |
| | Signal Bias (OSB).

Instructions on how to use phase bias corrections
                    provided for the final products:
					http://ftp.aiub.unibe.ch/CODE/IAR_README.TXT





### Meaning of BIA files
* Bias Solution INdependent EXchange Format (Bias-SINEX)
*-------------------------------------------------------------------------------
* CODE'S MGEX IAR phase/code OSB results for day 135, 2023     21-May-2023 02:10




### Center for Orbit Determination in Europe (CODE) and German Aerospace Center (DLR
CODE, the Center for Orbit Determination in Europe, is a joint venture of the following four institutions: Astronomical Institute, University of Bern (AIUB), Bern, Switzerland; Federal Office of Topography swisstopo, Wabern, Switzerland; Federal Agency of Cartography and Geodesy (BKG), Frankfurt a. M., Germany; Institut für Astronomische und Physikalische Geodäsie, Technische Universität München (IAPG, TUM), Munich, Germany.
It acts as a global analysis center of the International GNSS Service (IGS). The operational computations are performed at AIUB using the latest development version of the Bernese GNSS Software.
In this context a multi-GNSS solution is generated considering all active GPS, GLONASS, Galileo, BeiDou (expect for GEOs), and QZSS satellites as a contribution to the IGS-MGEX project. The results are published with a delay of about two weeks.


####### CLK Files
http://navigation-office.esa.int/products/gnss-products/
http://navigation-office.esa.int/products/gnss-products/esm.acn
from ESA: ESOC USING NAPEOS                                           ANALYSIS CENTER     
This document summarizes the Processing Strategy of the our published ESOC MGNSS products
Software used: NAPEOS 4.8
model description on http://navigation-office.esa.int/products/gnss-products/esm.pdf
Frequencies used:
GAL: L1-L5Q 
GPS: L1W-L2W (GPS C1/P2' code biases corrected to P1/P2 in GnssObs)

-> this means that CLK files contain iono free E1-E5a for GAL and L1W-L2W/L1P-L2P for GPS 




####### CODE software (CODE.ACN)
http://ftp.aiub.unibe.ch/CODE/CODE.ACN
INTERNATIONAL GNSS SERVICE     CODE Analysis Strategy Summary

 CLK: Clock corrections are consistent with          |
|                   |        carrier phase as well as P1/P2 pseudorange     |
|                   |        measurements.                                  |
|                   |        CODE P1-C1 pseudorange bias values of a        |
|                   |        moving 30-day solution are considered to       |
|                   |        correct C1/X2 and C1/P2 receiver data.

CLK Header
C  GPSEST V5.5        CODE.OSB @ ftp.aiub.unibe.ch/CODE/         SYS / DCBS APPLIED  
E  GPSEST V5.5        CODE.OSB @ ftp.aiub.unibe.ch/CODE/         SYS / DCBS APPLIED  
G  GPSEST V5.5        CODE.OSB @ ftp.aiub.unibe.ch/CODE/         SYS / DCBS APPLIED  
J  GPSEST V5.5        CODE.OSB @ ftp.aiub.unibe.ch/CODE/         SYS / DCBS APPLIED  
R  GPSEST V5.5        CODE.OSB @ ftp.aiub.unibe.ch/CODE/         SYS / DCBS APPLIED

so for GPS:
	CLK data is consistent with 1W/2W code and phase data. 
	The raw observables from the receiver 1C and 2C are transformed to 1W and 2W
	
The OSB file has DCBS for:
	GPS: C1C, C1L, C1X, C1W, C2L, C2S, C2X, C2W
    GAL: C1C, C1X, C5Q, C5X
	
	
	CLK: Clock corrections are consistent with carrier  |
|                   |        phase as well as P1/P2 pseudorange             |
|                   |        measurements.                                  |
|                   |        CODE P1-C1 pseudorange bias values of a moving |
|                   |        30-day solution are considered to correct      |
|                   |        C1/X2 and C1/P2 receiver data. 