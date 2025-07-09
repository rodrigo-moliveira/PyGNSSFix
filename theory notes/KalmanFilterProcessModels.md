# 📡 PPP EKF State Modeling and Noise Settings

This document describes the noise modeling assumptions and typical parameter values for an Extended Kalman Filter (EKF) implementation in a Precise Point Positioning (PPP) algorithm.

## 🧾 State Vector

The EKF considered in this project includes the following states:

- `position`
- `clock_bias`
- `ISB` (Inter-System Bias)
- `iono` (Ionosphere delay)
- `tropo` (Troposphere wet delay)
- `ambiguity` (float or fixed)
- `phase_bias`

---

## 🎯 Process Models and Noise Parameters


| **State**   | **Model**    | **Φₖ**           | **Qₖ**          | **Recommended Parameters**                               | **Units of σ**             |
|-------------|--------------|------------------|-----------------|----------------------------------------------------------|----------------------------|
| Position    | Random Walk  | `1`              | `σ² ⋅ Δt`       | σ ≈ 0.001–0.01 m/√s (static) or 0.1–1.0 m/√s (kinematic) | meters / √s                |
| Clock Bias  | Random Walk  | `1`              | `σ² ⋅ Δt`       | σ ≈ 1–3 m/√s → 3.3–10 ns/√s                              | seconds / √s               |
| ISB         | Random Walk  | `1`              | `σ² ⋅ Δt`       | σ ≈ 0.01–0.1 m/√s                                        | meters / √s                |
| Ionosphere  | Random Walk  | `1`              | `σ² ⋅ Δt`       | σ ≈ 0.1–1.0 m/√s                                         | meters / √s                |
| Ionosphere  | Gauss-Markov | `α = exp(-Δt/τ)` | `σ² ⋅ (1 - α²)` | σ ≈ 2.0 m, τ ≈ 600–1800 s                                | meters                     |
| Troposphere | Random Walk  | `1`              | `σ² ⋅ Δt`       | σ ≈ 0.001–0.01 m/√s                                      | meters / √s                |
| Troposphere | Gauss-Markov | `α = exp(-Δt/τ)` | `σ² ⋅ (1 - α²)` | σ ≈ 0.05–0.20 m, τ ≈ 3600–10800 s                        | meters                     |
| Ambiguity   | Random Walk  | `1`              | `σ² ⋅ Δt`       | σ ≈ 0.1–10 cycles/√s (depends on ambiguity resolution)   | cycles / √s or meters / √s |
| Phase Bias  | Random Walk  | `1`              | `σ² ⋅ Δt`       | σ ≈ 0.001–0.01 m/√s                                      | meters / √s                |


---

## 📌 Detailed Notes per State

### 1. Position
- **Static user**: use random walk with small noise (e.g., 1 cm/√s).
- **Kinematic user**: use either random walk or a constant velocity model (e.g., include velocity states).

### 2. Clock Bias
- GPS L1 time offset drift is about 10 ns/s ≈ 3 m/s.
- Random walk modeling with `σ_c = 1–3 m/√s` is common.
- You can optionally model **clock drift** as a separate state.

### 3. ISB (Inter-System Bias)
- Slowly varying between constellations (e.g., GPS–Galileo).
- Modeled as either a small random walk or a constant bias with low noise.

### 4. Ionosphere
- Only needed if **not using ionosphere-free combination**.
- Estimate one per satellite. Use random walk or Gauss-Markov (`τ ≈ 1000–2000 s`).
- You may scale `Q` with satellite elevation angle to improve realism.

### 5. Troposphere
- Zenith Wet Delay is modeled as a Gauss-Markov process.
- Typical parameters: `σ = 0.01 m`, `τ = 1800 s` (30 minutes).
- One state per receiver, not per satellite.

### 6. Ambiguities
- Modeled as constant for continuous tracking.
- When **float**, use small noise (e.g., `σ = 0.01–1 m`).
- When **fixed**, set process noise to zero and apply integer fix.

### 7. Phase Bias (OSB)
- Usually constant; can be initialized from IGS/CODE products.
- Typically not time-evolving; use `Φ = 1`, `Q = 0`.

---

## 🔧 EKF Q Matrix Construction

Your `Qₖ` (process noise covariance) should be block-diagonal, with blocks like:

```
Qₖ = diag(
σ_p² ⋅ Δt ⋅ I, # position
σ_c² ⋅ Δt, # clock bias
σ_ISB² ⋅ Δt, # ISB
σ_iono² ⋅ Δt ⋅ I_N, # iono per sat (if modeled)
σ_tropo² ⋅ (1 - α²), # tropo (GM)
... # ambiguities, etc.
)

```




## ✅ Typical Sigma Values

| **State**         | **σ (Static User)** | **σ (Kinematic User)** |
|-------------------|---------------------|------------------------|
| Position          | `0.01 m/√s`         | `0.05–0.1 m/√s`        |
| Clock Bias        | `1–3 m/√s`          | `1–3 m/√s`             |
| ISB               | `0.01 m/√s`         | `0.02 m/√s`            |
| Ionosphere        | `0.02 m/√s`         | `0.05 m/√s`            |
| Troposphere       | `0.01 m`            | `0.01 m`               |
| Ambiguity (float) | `0.01–1 m`          | `0.05–1 m`             |
| Phase Bias        | `0` (usually fixed) | `0`                    |

---


# 🔄 Process Noise Models in EKF

This section explains commonly used state evolution models in GNSS EKFs and how to construct their corresponding **State Transition Matrix** (`Φₖ`) and **Process Noise Covariance** (`Qₖ`) for time step `Δt`.

---

## 📌 1. Random Constant Model

**Description**:  
- Assumes the state does **not change over time**.
- Used for static biases or fixed parameters (e.g., phase bias, fixed ambiguities).

**State Equation**:  
```
xₖ = xₖ₋₁
```

**STM**:  
```
Φₖ = 1
```

**Process Noise**:  
```
Qₖ = 0  (or a very small value for numerical stability)
```

---


## 📌 2. White Noise Model

**Description**:  
- Each new state value is a sample from a **zero-mean Gaussian distribution**.
- No temporal correlation between states.
- Rarely used for states that evolve in time (too noisy), but can model sudden noise bursts.

**State Equation**:  
```
xₖ = wₖ
```

**STM**:  
```
Φₖ = 0
```

**Process Noise (already in discrete time)**:  
```
Qₖ = σ²
```

**Discretization Notes**:
- If you're using a **discrete-time white noise model directly**, do **not** multiply by `Δt`. `σ²` above is in units `unit²`
- If you're deriving from a **continuous-time noise power spectral density** `q` (in units of variance/sec), then discretize as:

```
Qₖ = q ⋅ Δt
```

This distinction is important for tuning based on physical noise characteristics.

**Units and Physical Meaning**:

Continuous-time white noise is defined via a **power spectral density (PSD)** `q`, with units:

```
[unit² / s]
```

When discretized over time step `Δt`, the process noise becomes:

```
Qₖ = q ⋅ Δt
```

and has units:

```
[unit²]
```

If you refer to the **square root of the PSD** (noise amplitude: square root of `q`), it has units:

```
[unit / √s]
```

This matches the **rate parameter** used in a **random walk** process noise:

`σ_RW` has units `[unit / √s]`

These unit relationships are important when tuning noise levels from physical models or sensor specs.
---

## 📌 3. Random Walk Model

**Description**:  
- Each state value is the sum of the previous state and a random perturbation.
- Used when the state changes slowly over time (e.g., position, clock bias, ionosphere delay).

**State Equation**:  
```
xₖ = xₖ₋₁ + wₖ
```

**STM**:  
```
Φₖ = 1
```

**Process Noise**:  
```
Qₖ = σ² ⋅ Δt
```

Where:
- `σ` (continuous-time sigma) is the **rate of change** (in units/√s).
- `Δt` is the time step (seconds).
- `wₖ` has units [unit]
- `Qₖ` has units [unit²]

**Multidimensional Extension** (e.g., 3D position):  
```
Φₖ = Iₙ  
Qₖ = σ² ⋅ Δt ⋅ Iₙ
```

---

## 📌 4. Gauss-Markov Model (1st Order)

**Description**:  
- Captures exponential decay in state value over time.
- Suitable for slowly varying quantities with mean-reverting behavior (e.g., troposphere, ionosphere).

**State Equation**:  
```
xₖ = α ⋅ xₖ₋₁ + wₖ
```

Where:
- `α = exp(-Δt / τ)`  
- `τ` is the **correlation time constant** (in seconds).

**STM**:  
```
Φₖ = α
```

**Process Noise**:  
```
Qₖ = σ² ⋅ (1 - α²)
```

Where `σ` is the **steady-state standard deviation**.

---

## 🧮 Summary Table

| **Model**                     | **Equation**         | **Φₖ**           | **Qₖ**          | Units of Sigma               |
|-------------------------------|----------------------|------------------|-----------------|------------------------------|
| Random Constant               | `xₖ = xₖ₋₁`          | `1`              | `0` or small ε  | -                            |
| White Noise (continuous-time) | `xₖ = wₖ`            | `0`              | `q ⋅ Δt`        | `q` has units `[unit² / s]`  |
| Random Walk                   | `xₖ = xₖ₋₁ + wₖ`     | `1`              | `σ² ⋅ Δt`       | `σ²` has units `[unit² / s]` |
| Gauss-Markov                  | `xₖ = α ⋅ xₖ₋₁ + wₖ` | `α = exp(-Δt/τ)` | `σ² ⋅ (1 - α²)` | `σ²` has units `[unit²]`     |

---

## ✅ Guidelines

- Use **Random Walk** for time-varying parameters like **position**, **clock bias**, and **ionosphere delay**.
- Use **Gauss-Markov** for slowly drifting but mean-reverting quantities like **troposphere delay**.
- Use **Random Constant** for fixed biases like **ambiguities (after fixing)** or **phase biases**.
- Use **White Noise** only for modeling bursts or if the state resets every time step (rare in PPP/INS/GNSS).

