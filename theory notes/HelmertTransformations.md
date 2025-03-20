# Helmert Transformation and ITRF Realization Transformations

## 1. Introduction
The **Helmert transformation** (also known as the **7-parameter transformation**) is used to transform coordinates between different reference frames. It consists of three **translation** parameters, three **rotation** parameters, and a **scale factor**. In the context of **International Terrestrial Reference Frame (ITRF) realizations**, this transformation is crucial for maintaining consistency between different ITRF versions (e.g., ITRF93 to ITRF2020).

## 2. Helmert Transformation Equations

The transformation from a source frame \( X_{\text{ITRF93}} \) to a target frame \( X_{\text{ITRF2020}} \) is given by:

\[
X_{\text{ITRF2020}} = T + (1 + S) R X_{\text{ITRF93}}
\]

where:
- \( X_{\text{ITRF93}} \) is the position vector in ITRF93 (meters)
- \( X_{\text{ITRF2020}} \) is the transformed position vector in ITRF2020 (meters)
- \( T = [T_x, T_y, T_z]^T \) is the translation vector (meters)
- \( S \) is the scale difference (unitless, parts per billion converted to a scale factor)
- \( R \) is the rotation matrix, which for small angles is approximated as:
  
  \[
  R \approx \begin{bmatrix} 1 & -R_z & R_y \\ R_z & 1 & -R_x \\ -R_y & R_x & 1 \end{bmatrix}
  \]

where \( R_x, R_y, R_z \) are the rotation parameters in radians.

## 3. Time-Dependent Transformation
Since reference frames evolve over time, the transformation parameters also change. The time-dependent transformation is given by:

\[
T(t) = T_{\text{ref}} + \dot{T} (t - t_{\text{ref}})
\]

\[
R(t) = R_{\text{ref}} + \dot{R} (t - t_{\text{ref}})
\]

\[
S(t) = S_{\text{ref}} + \dot{S} (t - t_{\text{ref}})
\]

where:
- \( t \) is the target epoch (e.g., 2025.0)
- \( t_{\text{ref}} \) is the reference epoch (e.g., 2015.0 for ITRF2020 transformations)
- \( \dot{T} = [\dot{T}_x, \dot{T}_y, \dot{T}_z]^T \) is the rate of change of translation (meters/year)
- \( \dot{R} = [\dot{R}_x, \dot{R}_y, \dot{R}_z]^T \) is the rate of change of rotation (radians/year)
- \( \dot{S} \) is the rate of change of scale difference (unitless/year)

These equations ensure that the transformation accounts for the time evolution of reference frames.

## 4. Application to ITRF93 to ITRF2020 Transformation
Using the transformation parameters provided by the **EPSG:9998** transformation model:

**Reference parameters (at epoch 2015.0):**
- \( T_{\text{ref}} = [65.8, -1.9, 71.3] \) mm
- \( R_{\text{ref}} = [3.36, 4.33, -0.75] \) mas
- \( S_{\text{ref}} = -4.47 \) ppb

**Rates of change:**
- \( \dot{T} = [2.8, 0.2, 2.3] \) mm/year
- \( \dot{R} = [0.11, 0.19, -0.07] \) mas/year
- \( \dot{S} = -0.12 \) ppb/year

### **Unit Conversions**
- **Translation:** 1 mm = **0.001 meters**
- **Rotation:** 1 mas = **4.84813681109536e-09 radians**
- **Scale difference:** 1 ppb = **1e-9 (unitless)**

### **Example Computation for Epoch 2025.0**
Computing the adjusted parameters at **t = 2025.0**:

\[
T(2025) = T_{\text{ref}} + \dot{T} (2025 - 2015)
\]
\[
R(2025) = R_{\text{ref}} + \dot{R} (2025 - 2015)
\]
\[
S(2025) = S_{\text{ref}} + \dot{S} (2025 - 2015)
\]

Once these adjusted parameters are computed, the Helmert transformation is applied to the input position vector.

## 5. Conclusion
The **Helmert transformation** is a critical tool for converting coordinates between ITRF realizations. Given the ongoing updates in terrestrial reference frames, the time-dependent transformation equations ensure **accurate positioning** over time.

By applying these transformations, consistency between datasets using different ITRF versions (such as CSpice with ITRF93 and modern GNSS products using ITRF2020) can be maintained.

For more details on the Helmert transformation and ITRF realizations, refer to the official IERS documentation and transformation models like EPSG:9998.