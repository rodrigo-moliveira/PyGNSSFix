# EKF Variants with Pseudorange and Range Rate

## State Vector

\[
\mathbf{x} =
\begin{bmatrix}
\mathbf{p} \\
\mathbf{v} \\
c \cdot dt \\
c \cdot \dot{dt}
\end{bmatrix}
\]

---

## Case 1: EKF with PV Model (PR only)

**Prediction Step**:

\[
F =
\begin{bmatrix}
I & \Delta t \cdot I & 0 & 0 \\
0 & I & 0 & 0 \\
0 & 0 & 1 & \Delta t \\
0 & 0 & 0 & 1
\end{bmatrix}
,\quad
\mathbf{x}_{k|k-1} = F \cdot \mathbf{x}_{k-1|k-1}
\]

**Update Step**:

Pseudorange measurement model:

\[
z_{\text{PR}} = \|\mathbf{p}^{\text{sat}} - \mathbf{p}\| + c \cdot dt
\]

Jacobian \( H \):

\[
H =
\begin{bmatrix}
\frac{\partial \rho}{\partial \mathbf{p}} & 0 & 1 & 0
\end{bmatrix}
\]

---

## Case 2: EKF with PV Model (PR + Range Rate)

**Prediction Step**:

Same as above.

**Update Step**:

Pseudorange: same as Case 1.

Range Rate measurement:

\[
z_{\text{RR}} = \left( \mathbf{v}^{\text{sat}} - \mathbf{v} \right) \cdot \hat{\boldsymbol{\rho}} + c \cdot \left( \dot{dt}^{\text{rec}} - \dot{dt}^{\text{sat}} - \text{rel}_{\text{sat}} \right)
\]

Jacobian \( H \):

\[
H =
\begin{bmatrix}
\frac{\partial \rho}{\partial \mathbf{p}} & 0 & 1 & 0 \\
0 & -\hat{\boldsymbol{\rho}}^T & 0 & 1
\end{bmatrix}
\]

Where \( \hat{\boldsymbol{\rho}} = \frac{\mathbf{p}^{\text{sat}} - \mathbf{p}}{\|\mathbf{p}^{\text{sat}} - \mathbf{p}\|} \) is the line-of-sight unit vector.

---

## Case 3: EKF with Random Walk (PR + Range Rate)

**Prediction Step**:

\[
F = I
,\quad
\mathbf{x}_{k|k-1} = \mathbf{x}_{k-1|k-1}
\]

**Update Step**:

Use same pseudorange and range rate equations as in Case 2.

---

## Clock Bias and Drift Model

State:

\[
\begin{bmatrix}
c \cdot dt \\
c \cdot \dot{dt}
\end{bmatrix}
\]

**Prediction**:

\[
\begin{bmatrix}
c \cdot dt \\
c \cdot \dot{dt}
\end{bmatrix}_{k|k-1}
=
\begin{bmatrix}
1 & \Delta t \\
0 & 1
\end{bmatrix}
\cdot
\begin{bmatrix}
c \cdot dt \\
c \cdot \dot{dt}
\end{bmatrix}_{k-1|k-1}
\]

**Process Noise Covariance**:

\[
Q_{\text{clock}} = \sigma^2
\begin{bmatrix}
\frac{\Delta t^3}{3} & \frac{\Delta t^2}{2} \\
\frac{\Delta t^2}{2} & \Delta t
\end{bmatrix}
\]