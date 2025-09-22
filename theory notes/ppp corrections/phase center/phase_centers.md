# GNSS Receiver Phase Center Offsets (PCO) and Phase Center Variations (PCV) in PPP

The distance between a GNSS satellite and a receiver is measured between their
respective Antenna Phase Centers (APCs). If the precise satellite orbit data are referred to
the satellite Center of Mass (CoM), such as in the IGS SP3 files and some real-time products,
the difference between the satellite CoM and APC must be accounted for. On the other
hand, the receiver’s reference location is often determined at a marker position that can
be the Antenna Reference Point (ARP) or another point on the ground below the receiver
antenna. Consequently, the difference between the receiver reference point and its APC
must be considered.
APC errors are frequency-dependent and can be classified into Phase Center Offset (PCO) and Phase Center Variation (PCV). Calibrated corrections for these errors are
available in the ANTenna Exchange (ANTEX) format file provided by the IGS

## 1. Introduction

GNSS measurements are referenced to the so-called antenna phase center, which exists for both the satellite and the receiver. The phase center is not necessarily the geometric center of the antenna, nor is it constant. Instead, it depends on the direction of the incoming radio signal.

The phase center is defined as the apparent source of radiation. In an ideal case, it would have a spherical equiphase contour, but real-world antennas exhibit irregularities, causing variations in the apparent radiation origin.

Furthermore, the phase center of an antenna is dependent on azimuth, elevation, and frequency. A simple model assumes that phase centers differ only along the vertical axis of the antenna.

Precise Point Positioning (PPP) requires accurate modeling of GNSS measurements, including corrections for receiver antenna phase center offsets (PCO) and phase center variations (PCV). This document details the methodology for incorporating these corrections into GNSS PPP algorithms.

## 2. Reference Points and Coordinate Systems

For a GNSS receiver antenna, the following reference points are defined:

- **Marker Position (r_marker):** The geodetic reference monument on which an antenna is mounted directly with forced centering or on a tripod, typically reported in IGS solutions.
- **Antenna Reference Point (ARP) (r_ARP):** A well-defined point on the antenna, such as the center of the bottom surface of the preamplifier. The antenna height is measured from the marker to the ARP and reported in the ANTENNA: DELTA H/E/N header record. Small horizontal eccentricities of the ARP relative to the marker can also be reported in this record.
- **Antenna Phase Center (APC) (r_APC):** The effective reception point of the signal, which depends on frequency, azimuth, and elevation. It is frequency-dependent and varies with the minimum elevation angle.

The relationships between these positions are:

```
r_ARP = r_marker + Delta_ARP
r_APC = r_ARP + PCO(nu)
```

where:

- **Delta_ARP** is the offset from the marker to the ARP (provided in ENU coordinates in ANTEX files),
- **PCO(nu)** is the frequency-dependent phase center offset (also in ENU coordinates),
- **nu** represents the frequency of the GNSS signal.

## 3. Receiver Antenna Phase Center

GNSS measurements are referenced to the antenna phase center (APC). Since this location is frequency-dependent, the antenna reference point (ARP) is used as a more stable reference. Manufacturers provide technical information on the APC position relative to the ARP, and IGS compiles APC corrections in PCV and ANTEX files for various antenna models.



### PCO Correction

The phase center offset (PCO) is a position vector that must be projected onto the receiver-to-satellite line-of-sight (LOS) vector **e** to obtain a range correction:

```
zeta_PCO,r,j = - e ⋅ (A r_PCO,r,j)
```

where:

- **e** is the unit LOS vector from receiver to satellite,
- **r_PCO,r,j** is the PCO vector in the antenna-fixed coordinate system,
- **A** is the direction cosine matrix that transforms antenna-fixed coordinates to ECEF coordinates,
- The negative sign ensures the correction is applied in the proper direction.

The corrected pseudorange and carrier phase measurements include this correction:

```
P_tilde = P + zeta_PCO,r,j
Phi_tilde = Phi + zeta_PCO,r,j
```

### PCV Correction

Phase center variations (PCV) represent the variation of the phase center as a function of azimuth and elevation and are provided in ANTEX files as scalar values in meters. PCVs are directly applied as range corrections:

```
P_tilde = P + zeta_PCO,r,j + PCV(alpha, theta, nu)
Phi_tilde = Phi + zeta_PCO,r,j + PCV(alpha, theta, nu)
```

where:

- **PCV(alpha, theta, nu)** is interpolated from ANTEX data using the satellite azimuth **alpha** and elevation **theta**.

## 4. Extracting Data from OBS and ANTEX Files

### From Observation (OBS) Files:

- **APPROX POSITION XYZ** → Geocentric approximate marker position (**r_M**)
- **ANTENNA: DELTA H/E/N** → Antenna height and eccentricity of ARP relative to marker (**delta_ARP**)

### From ANTEX Files:

- **NORTH / EAST / UP:** Receiver antenna eccentricities of the mean antenna phase center relative to the ARP.
- **Non-azimuth-dependent pattern (NOAZI):** Phase pattern values in millimeters from ZEN1 to ZEN2.
- **Azimuth-dependent pattern (DAZI > 0.0):** Phase pattern values provided for each azimuth angle, followed by values from ZEN1 to ZEN2.


## 5. References

- IGS ANTEX Format Documentation. https://files.igs.org/pub/data/format/antex14.txt
- Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck, Springer Cham, 2017
- ESA GNSS DATA PROCESSING, Volume I: Fundamentals and Algorithms, J. Sanz Subirana, J.M. Juan Zornoza and M. Hernández-Pajares


# GNSS Receiver Phase Center Offsets (PCO) and Phase Center Variations (PCV) in PPP

## 1. Introduction

GNSS measurements are referenced to the so-called antenna phase center, which exists for both the satellite and the receiver. The phase center is not necessarily the geometric center of the antenna, nor is it constant. Instead, it depends on the direction of the incoming radio signal.

The phase center is defined as the apparent source of radiation. In an ideal case, it would have a spherical equiphase contour, but real-world antennas exhibit irregularities, causing variations in the apparent radiation origin.

Furthermore, the phase center of an antenna is dependent on azimuth, elevation, and frequency. A simple model assumes that phase centers differ only along the vertical axis of the antenna.

Precise Point Positioning (PPP) requires accurate modeling of GNSS measurements, including corrections for receiver antenna phase center offsets (PCO) and phase center variations (PCV). This document details the methodology for incorporating these corrections into GNSS PPP algorithms.

## 2. Reference Points and Coordinate Systems

For a GNSS receiver antenna, the following reference points are defined:

- **Marker Position (r_marker):** The geodetic reference monument on which an antenna is mounted directly with forced centering or on a tripod, typically reported in IGS solutions.
- **Antenna Reference Point (ARP) (r_ARP):** A well-defined point on the antenna, such as the center of the bottom surface of the preamplifier. The antenna height is measured from the marker to the ARP and reported in the ANTENNA: DELTA H/E/N header record. Small horizontal eccentricities of the ARP relative to the marker can also be reported in this record.
- **Antenna Phase Center (APC) (r_APC):** The effective reception point of the signal, which depends on frequency, azimuth, and elevation. It is frequency-dependent and varies with the minimum elevation angle.

The relationships between these positions are:

```
r_ARP = r_marker + Delta_ARP
r_APC = r_ARP + PCO(nu)
```

where:

- **Delta_ARP** is the offset from the marker to the ARP (provided in ENU coordinates in ANTEX files),
- **PCO(nu)** is the frequency-dependent phase center offset (also in ENU coordinates),
- **nu** represents the frequency of the GNSS signal.

### Receiver Antenna Phase Center and Antenna Reference Point

GNSS measurements are referenced to the antenna phase center (APC). Since this location is frequency-dependent, the antenna reference point (ARP) is used as a more stable reference. Manufacturers provide technical information on the APC position relative to the ARP, and IGS compiles APC corrections in PCV and ANTEX files for various antenna models.

In geodetic positioning, receiver coordinates are referenced to a monument marker (MM) or an external benchmark (BM). IGS SINEX files provide:

- **MM coordinates** in ECEF coordinates (SOLUTION/ESTIMATE block),
- **ARP relative to MM** in UNE coordinates (SITE/ECCENTRICITY block),
- **APC offsets** for different frequencies and calibration models (SITE/GPS PHASE CENTER block).

## 3. PCO Correction

The phase center offset (PCO) is a position vector that must be projected onto the receiver-to-satellite line-of-sight (LOS) vector **e** to obtain a range correction:

```
zeta_PCO,r,j = - e ⋅ (A r_PCO,r,j)
```

where:

- **e** is the unit LOS vector from receiver to satellite,
- **r_PCO,r,j** is the PCO vector in the antenna-fixed coordinate system,
- **A** is the direction cosine matrix that transforms antenna-fixed coordinates to ECEF coordinates,
- The negative sign ensures the correction is applied in the proper direction.

The corrected pseudorange and carrier phase measurements include this correction:

```
P_tilde = P + zeta_PCO,r,j
Phi_tilde = Phi + zeta_PCO,r,j
```

## 4. PCV Correction

Phase center variations (PCV) represent the variation of the phase center as a function of azimuth and elevation and are provided in ANTEX files as scalar values in meters. PCVs are directly applied as range corrections:

```
P_tilde = P + zeta_PCO,r,j + PCV(alpha, theta, nu)
Phi_tilde = Phi + zeta_PCO,r,j + PCV(alpha, theta, nu)
```

where:

- **PCV(alpha, theta, nu)** is interpolated from ANTEX data using the satellite azimuth **alpha** and elevation **theta**.

## 5. Extracting Data from OBS and ANTEX Files

### From Observation (OBS) Files:
- **APPROX POSITION XYZ** → Geocentric approximate marker position (**r_M**)
- **ANTENNA: DELTA H/E/N** → Antenna height and eccentricity of ARP relative to marker (**delta_ARP**)

### From ANTEX Files:
- **NORTH / EAST / UP:** Receiver antenna eccentricities of the mean antenna phase center relative to the ARP.
- **Non-azimuth-dependent pattern (NOAZI):** Phase pattern values in millimeters from ZEN1 to ZEN2.
- **Azimuth-dependent pattern (DAZI > 0.0):** Phase pattern values provided for each azimuth angle, followed by values from ZEN1 to ZEN2.

## 6. Implementation in a PPP Algorithm

1. **State Representation:** The estimated receiver state should be the marker position **r_marker**.
2. **Observation Reconstruction:**
   - Compute **r_ARP** and **r_APC** from known offsets.
   - Correct raw pseudorange and carrier phase measurements using PCO.
   - Apply PCV correction as an additional range correction.
3. **Handling in Observables Domain:**
   - **PCO is a position offset**, so it must be **projected onto the LOS vector**.
   - **PCV is a range correction**, so it is **directly added**.

By applying these corrections in the observation domain, the estimated marker position remains consistent, ensuring accurate PPP processing.

## 7. References

- IGS ANTEX Format Documentation
- GNSS Textbooks on PPP and Antenna Corrections

## 8. Choosing the Reference Position for Range Computation

To ensure consistency in PPP estimation, the **marker position** (r_marker) should be used to compute the geometric range and the line-of-sight (LOS) vector:

1. **Geometric range calculation**:
   - Compute the range between **r_marker** and the satellite position.
   - This ensures that the estimated position in PPP refers directly to the marker.

2. **LOS vector computation**:
   - Compute the unit vector from **r_marker** to the satellite.
   - The same LOS vector is then used for all frequency-dependent PCO corrections.

3. **Impact on Estimation**:
   - Using **r_marker** ensures that the estimated state corresponds to the geodetic marker.
   - If **r_ARP** were used instead, the estimated position would correspond to **r_ARP**, requiring post-processing corrections.

Thus, using the marker position from the start simplifies the modeling and ensures that the final estimated position remains consistent with the reference geodetic monument.

## 9. 8. Reconstruction of Observables in PPP

To reconstruct the GNSS observables (pseudorange and carrier phase) while accounting for phase center offsets (PCO) and phase center variations (PCV), we follow a structured approach.
8.1 Geometric Range Computation

The geometric range ρ is computed using the marker position r_marker:
ρ = ||r_sat − r_marker||

Where:

    r_sat is the satellite position in ECEF,
    r_marker is the receiver marker position in ECEF.

The line-of-sight (LOS) unit vector is then:

e = (r_sat − r_marker) / ρ

Using the marker position ensures consistency across all signals and allows reusing the computed LOS vector for different frequencies.

8.2 Applying Antenna Corrections
8.2.1 Correction for Mounting Offset (Marker → ARP)

The offset from marker to antenna reference point (ARP), ΔARP, is provided in ENU coordinates in the RINEX observation file. It must be transformed to ECEF coordinates and added to the marker position:

r_ARP = r_marker + A_ENU→ECEF * ΔARP

Where A_ENU→ECEF is the transformation matrix from ENU to ECEF at the receiver's location.

8.2.2 Correction for Phase Center Offset (PCO)

The phase center offset (PCO) r_PCO is provided in ENU coordinates for each frequency. Like ΔARP, it must be transformed into ECEF and added to compute the phase center position:

r_APC = r_ARP + A_ENU→ECEF * r_PCO(f)

Where r_PCO(f) is frequency-dependent.

The projection of PCO onto the LOS vector gives the range correction:

ζ_PCO,r,j = −e ⋅ (A_ENU→ECEF * r_PCO,r,j)

8.2.3 Correction for Phase Center Variation (PCV)

The PCV is already in the range domain and is directly interpolated from the ANTEX file as a function of azimuth and elevation:

ζ_PCV,r,j = PCV(α,θ,f)

Where:

    α is the azimuth of the satellite as seen from the receiver,
    θ is the elevation angle,
    f is the signal frequency.

8.3 Corrected Observables

The corrected pseudorange and carrier phase observables are:

P' = P + ζ_PCO,r,j + ζ_PCV,r,j

Φ' = Φ + ζ_PCO,r,j + ζ_PCV,r,j

Where:

    P and Φ are the raw pseudorange and carrier phase measurements, respectively.
    ζ_PCO,r,j is the PCO correction (computed as a projection).
    ζ_PCV,r,j is the PCV correction (directly taken from ANTEX).

This approach ensures that all observations are referenced correctly to the antenna phase center, allowing for consistent and accurate PPP processing.


# Satellite APC Errors
The PCO of a satellite is defined as the distance between its CoM and the mean APC,
expressed in the satellite body-fixed frame. The IGS definition for this frame is that the
coordinate system’s origin is at the satellite’s CoM, the z-axis is parallel to the antenna
boresight, the y-axis is aligned with the rotation axis of the solar panels, and the x-axis
completes the right-handed system