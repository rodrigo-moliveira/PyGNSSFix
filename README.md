# PyGNSSFix

**PyGNSSFix** is a Python-based GNSS positioning toolkit implementing several algorithms for post-processing analysis.  
Currently, the available modes include:

- **Single Point Positioning (SPP)**
- **Pseudorange-based Precise Point Positioning (PR-PPP)**
- **Carrier Phase-based Precise Point Positioning (CP-PPP)**  
  (with cycle slip detection and ambiguity resolution available)

Currently supported constellations: **GPS** and **Galileo**.


## 1. Features
- Flexible support for multiple GNSS constellations (GPS, GAL).
- RINEX-based processing workflow (post-processing positioning solution).
- Implements standard and precise GNSS positioning models.
- Post-processing scripts for result analysis, residuals and visualization.
- Highly configurable and extensible system.
- Two solvers are implemented:
  - Weighted Least Squares
  - Extended Kalman Filter


## 2. Algorithm Details \& Variations

**PyGNSSFix** supports several variations of the positioning algorithms (SPP, PR-PPP, CP-PPP) depending on:

1. **Constellations**  
2. **Observables**  
3. **Estimation type (Position / Velocity)**  
4. **Filtering method**  
5. **Stochastic process models**


### 2.1 Constellations

| Mode   | Supported Constellations |
|--------|--------------------------|
| SPP    | GPS, GAL, or GPS+GAL     |
| PR-PPP | GPS, GAL, or GPS+GAL     |
| CP-PPP | GPS, GAL, or GPS+GAL     |


### 2.2 Observables

| Observation Type         | Description                                                                               |
|--------------------------|-------------------------------------------------------------------------------------------|
| **Raw Single-Frequency** | Use pseudorange and/or carrier phase from a single frequency per constellation            |
| **Raw Dual-Frequency**   | Use pseudorange and/or carrier phase from two frequencies per constellation               |
| **Iono-Free**            | Linear combination of dual-frequency observations to remove first-order ionospheric delay |


### 2.3 Estimation Variables

- **Position-only**:  
  Estimate receiver position and clock bias (using PR/CP observations).  

- **Position + Velocity**:  
  If Doppler observations are enabled, estimate velocity and clock bias rate in addition to position and clock bias.


The process varies depending on the selected solver (WLS or EKF).

#### 2.3.1 Weighted Least Squares (WLS)

Two-step process:
1. First solve **position + clock bias** from PR/CP observations (Iterated Least Squares).  
2. Then solve **velocity + clock bias rate** from Doppler observations (separate LS problem).  


#### 2.3.2 Extended Kalman Filter (EKF)

All states are estimated **jointly** in the EKF state vector:

- Position  
- Velocity  
- Clock bias  
- Clock bias rate  

**State transitions** can be configured as:
- **Uncorrelated**: In this case the user can select different stochastic processes (see subsection below).  
- **Correlated**: For Position / Velocity and Clock Bias / Clock Bias Rate it is possible to correlate the states via dynamical models (PV models, clock models).  

For more information about the states, see the specific internal documentation [KalmanFilterProcessModels.md](theory%20notes/KalmanFilterProcessModels.md).


### 2.4 Stochastic Process Models

Each state variable (position, velocity, clock bias, clock bias rate) can evolve according to:

| Process Type     | Description                                         |
|------------------|-----------------------------------------------------|
| **White Noise**  | State considered uncorrelated between epochs        |
| **Random Walk**  | State evolves as cumulative noise over time         |
| **Gauss–Markov** | First-order autoregressive model with finite memory |

For more information about the stochastic processes, see the specific internal [KalmanFilterProcessModels.md](theory%20notes/KalmanFilterProcessModels.md).

### 2.5 Estimated States

The list of all available states that are estimable are the following:

* **Position** 
* **Velocity** 
* **Clock Bias**
* **Clock Bias Rate**
* **Troposphere Wet Delay**
* **Inter System Bias** 
* **Ionosphere Delay** 
* **Ambiguity**
* **Receiver Phase Bias**

More details in the table below:

| **State Name**        | **Dimension (with description)**                                   | **Unit** | **State description (availability)**                                               |
|-----------------------|--------------------------------------------------------------------|----------|------------------------------------------------------------------------------------|
| Position Vector       | 3×1 (ECEF)                                                         | m        | Mandatory state, always included                                                   |
| Velocity              | 3×1 (ECEF)                                                         | m/s      | Enabled by configuration, requires Doppler measurements                            |
| Clock Bias            | 1×1                                                                | s        | Mandatory state, always included                                                   |
| Clock Bias Rate       | 1×N_const (per constellation)                                      | s/s      | Enabled by configuration, requires Doppler measurements                            |
| Troposphere Wet Delay | 1×1                                                                | m        | Optional, included when troposphere modeling is enabled                            |
| Inter-System Bias     | 1×1                                                                | s        | Only estimated when dual constellations are processed (e.g., GPS + Galileo)        |
| Ionosphere Delay      | N_sat×1 (per satellite)                                            | m        | Estimated when raw dual-frequency observations are processed and enabled by config |
| Ambiguity             | N_obs×1 (per carrier-phase observation, per satellite & frequency) | cycles   | Only present in CP-PPP mode                                                        |
| Receiver Phase Bias   | N_freq×1 (per frequency, per constellation)                        | s        | Only present in CP-PPP mode                                                        |

### 2.6 Observation Models

The implemented measurement models are:

- **Pseudorange (PR)**
- **Carrier Phase (CP)**
- **Pseudorange Rate (PRR)**


**Single-Frequency Pseudorange (PR):**

$$
PR_j = \rho(t, t-\tau) + c \, \Delta t_{rec} + ISB - c \, (\Delta t_{sat} - b_{sat,j}) + \mu_j \, (I + \delta I) + T + pcc_{rec,j} + pcc_{sat,j}
$$

**Single-Frequency Carrier Phase (CP):**

$$
CP_j = \rho(t, t-\tau) + c \, \Delta t_{rec} + b^{\phi}_j + ISB - c \, (\Delta t_{sat} - b_{sat,j}) + \mu_j \, (I + \delta I) + T + \lambda_j N_j + pcc_{rec,j} + pcc_{sat,j}
$$

**Iono-Free Pseudorange (IF-PR):**

$$
PR_{IF} = \rho(t, t-\tau) + c \, \Delta t_{rec} + ISB - c \, (\Delta t_{sat} - b_{sat,IF}) + T + pcc_{rec,IF} + pcc_{sat,IF}
$$

**Iono-Free Carrier Phase (IF-CP):**  

$$
CP_{IF} = \rho(t, t-\tau) + c \, \Delta t_{rec} + b^{\phi}_{IF} + ISB - c \, (\Delta t_{sat} - b_{sat,IF}) + T + \lambda_{IF} N_{IF} + pcc_{rec,IF} + pcc_{sat,IF}
$$

**Pseudorange Rate (PRR):**  

$$
PRR_j = ( \mathbf{v}^{sat} - \mathbf{v}_{rec} ) \cdot \mathbf{los} + c \, (\dot{\Delta t}_{rec} - \dot{\Delta t}_{sat} - \dot{rel}_{sat})
$$


#### 2.6.1 List of Variables

- $\rho(t, t-\tau)$: Geometric range between receiver and satellite at transmit time  
- $c$: Speed of light  
- $\Delta t_{rec}$: Receiver clock bias, with hardware bias included (s)  
- $\dot{\Delta t}_{rec}$: Receiver clock drift (s/s)  
- $\Delta t_{sat}$: Satellite clock bias, with hardware bias included (depending on configuration) (s)  
- $\dot{\Delta t}_{sat}$: Satellite clock drift (s/s)  
- $\dot{rel}_{sat}$: Relativistic satellite clock correction rate (s/s)  
- $ISB$: Inter-system bias (between GNSS constellations)  
- $b_{sat,j}$: Satellite group delay (code bias) for frequency $j$  
- $b^{\phi}_j$: Receiver phase bias for frequency $j$  
- $I$: Ionospheric delay (a-priori model)  
- $\delta I$: Estimated ionospheric delay (residual)  
- $\mu_j$: Frequency-dependent ionospheric scaling factor  
- $T$: Tropospheric delay  
- $N_j$: Carrier-phase ambiguity for frequency $j$ (float or integer)  
- $\lambda_j$: Carrier wavelength for frequency $j$  
- $pcc_{rec,j}$: Receiver phase center correction (frequency-dependent)  
- $pcc_{sat,j}$: Satellite phase center correction (frequency-dependent)  
- $PR_j$: Pseudorange measurement at frequency $j$  
- $CP_j$: Carrier-phase measurement at frequency $j$  
- $PR_{IF}$: Iono-free pseudorange combination  
- $CP_{IF}$: Iono-free carrier-phase combination  
- $\mathbf{v}^{sat}$: Satellite velocity vector (ECEF)  
- $\mathbf{v}_{rec}$: Receiver velocity vector (ECEF)  
- $\mathbf{los}$: Line-of-sight unit vector from receiver to satellite  
- $PRR_j$: Pseudorange-rate (Doppler) observable at frequency $j$  


See [1] (Chapter 21) for a complete explanation of the meaning of all estimation variables and their relation with the 
theoretical terms that form the standard definition of the GNSS observation equations. 

For more details see the theory notes.

### 2.7 Corrections \& Models

The following correction models are implemented:

- **Ionospheric Delay**:  
  - Klobuchar Model
  - NTCM-G Model
  - Global Ionospheric Maps (IONEX)  

- **Tropospheric Delay**:  
  - Saastamoinen model  
  - GPT3 

- **Satellite Clock Corrections**:  
  - Broadcast clock parameters  
  - Precise clock files (CLK)  

- **Satellite Orbits**:  
  - Broadcast ephemerides  
  - Precise orbits (SP3)  

- **Antenna Phase Center Offsets/Variations**:  
  - Receiver antenna models from ANTEX  
  - Satellite antenna models from ANTEX 

For more details see the theory notes. 

## 3. Input Data

**PyGNSSFix** requires different sets of input files depending on the processing mode (SPP, PR-PPP, CP-PPP).  
The following file types are supported:


### 3.1 Core GNSS Files

- **Observation RINEX (`.obs`)**  
  Contains raw pseudorange, carrier-phase, doppler and C/N0 measurements collected by the receiver.  
  Used as the primary input for all positioning modes.

- **Navigation RINEX (`.nav`)**  
  Broadcast ephemerides, satellite clock corrections, and ionospheric parameters.  
  Required for Single Point Positioning (SPP).

- **Clock Files (`.clk`)**  
  Precise satellite clock corrections provided by analysis centers (e.g., IGS).  
  Essential for Precise Point Positioning (PPP).

- **Orbit Files (`.sp3`)**  
  Precise satellite orbit products (SP3c/d format).  
  Used in PR-PPP and CP-PPP for higher accuracy than broadcast orbits.

- **Ionosphere Maps (`.ionex`)**  
  Global ionospheric models in IONEX format, typically from IGS.  
  Provide vertical total electron content (VTEC) for ionospheric delay corrections.  

- **ANTEX Files (`.atx`)**  
  Antenna phase center offsets (PCO) and variations (PCV) for both satellites and receivers.  
  Required for sub-centimeter accuracy in PPP.

- **Bias Files (SINEX BIAS)**  
  Contain observable-specific biases (OSBs) or differential code biases (DCBs).  
  OSBs are required for ambiguity resolution in CP-PPP model (PPP-AR).


### 3.2 Auxiliary Data Files

- **CSpice Kernels**  
  External ephemeris and reference data used for high-precision transformations and for the position of the Sun and Moon.
  The required files are already contained in the repository folder ([cspice_kernels](workspace/geo_time_data/cspice_kernels))
  and only need to be changed if an update of the kernels is intended.
  - [de421.bsp](workspace/geo_time_data/cspice_kernels/de421.bsp) – Solar system planetary ephemerides  
  - [earth_latest_high_prec.bpc](workspace/geo_time_data/cspice_kernels/earth_latest_high_prec.bpc) – Earth orientation and rotation model  
  - [naif0012.tls.pc](workspace/geo_time_data/cspice_kernels/naif0012.tls.pc) – Leap seconds and time system definitions

- **Fault Injection Files**  
  Custom input defining artificial errors or anomalies to be injected into the measurements or models.  
  Useful for testing robustness and fault detection algorithms.
  An example is provided in [test_fault_config.csv](workspace/faults/test_fault_config.csv)

- **Frame Conversion Files (EPSG JSON)**  
  EPSG-based coordinate system definitions in JSON format.  
  Required for transformations between different terrestrial/space frames.
  The required files necessary for the current implementations of the software are already available in the repository folder (workspace/geo_time_data/ITRF)
  and only need to be changed if an update of the transformations is intended:
  - [Transformation: ITRF93 to ITRF2014](workspace/geo_time_data/ITRF/EPSG_8074.json)
  - [Transformation: ITRF93 to ITRF2020](workspace/geo_time_data/ITRF/EPSG_9998.json)

- **Earth Orientation Parameters (EOP)**  
  Files providing Earth orientation and time system corrections:  
  - [finals1980.all](workspace/geo_time_data/eop/finals1980.all) – IERS EOP series for the 1980 model-based transformations   
  - [finals2000.all](workspace/geo_time_data/eop/finals2000.all) – IERS EOP series for the 2000 model-based transformations
  - [tai-utc.dat](workspace/geo_time_data/eop/tai-utc.dat) – Leap seconds and UTC–TAI offset history

- **Tropospheric Model Grid Files (GPT)**  
  [gpt3_1.grd](workspace/geo_time_data/tropo/gpt3_5.grd) (Global Pressure and Temperature 3 model).  
  Provide site-dependent meteorological parameters for advanced tropospheric delay modeling.

### 3.3 Required Files per Mode

| **Mode** | **Mandatory Files**                                                                                       | **Optional Files**             |
|----------|-----------------------------------------------------------------------------------------------------------|--------------------------------|
| SPP      | Observation RINEX, Navigation RINEX                                                                       | –                              |
| PR-PPP   | Observation RINEX, SP3 (orbits), CLK (satellite clocks), SINEX Bias (with OSBs or DCBs)                   | Navigation RINEX, ANTEX, IONEX |
| CP-PPP   | Observation RINEX, SP3 (orbits), CLK (satellite clocks), SINEX Bias (**must be OSBs when AR is enabled**) | Navigation RINEX, ANTEX, IONEX |


## 4. Output files

The full list of output files available is:

* Output States:
  * `ambiguity.txt` : estimated ambiguity states
  * `receiver_phase_bias.txt` : estimated receiver phase biases states
  * `position.txt` : estimated position state
  * `velocity.txt` : estimated velocity state
  * `clock_bias.txt` : estimated clock bias
  * `iono.txt` : estimated ionospheric delay states
  * `isb.txt` : estimated ISB state
  * `tropo.txt` : estimated tropospheric wet delay state
  * `clock_bias_rate.txt` : estimated clock bias rate states
  
* Additional Files:
  * `satellite_azel.txt` : satellite azimuth and elevation angles
  * `DOP_ECEF.txt` : DOPs in ECEF frame
  * `DOP_ENU.txt` : DOPs in ENU frame
  * `prefit_residuals.txt` : estimation prefit residuals
  * `postfit_residuals.txt` : estimation postfit residuals
  * `time.txt` : time in string format
  * `observations.txt` : observations used in the estimation processed
  * `raw_observations.txt` : raw input observations
  * `melbourne_wubbena_obs.txt` : Melbourne Wubbena observations (combination)
  * `geometry_free_obs.txt` : Geometry-Free observations (combination)
  * `cycle_slips.txt` : detected cycle slips

Note that the files are only created when the data is processed by the selected configuration.


---

## 5. Installation

### 5.1 Requirements
- Python 3.10+
- Requirement libraries in [requirements.txt](requirements.txt). 

### 5.2 Setup

```bash
git clone https://github.com/rodrigo-moliveira/PyGNSSFix.git
cd PyGNSSFix
pip install -r requirements.txt
```

## 6. Usage

There are two main entry points:

### 6.1 Run the GNSS Algorithms

```bash
python algs/main_gnss.py <config_gnss.json>
```

* Executes the positioning algorithm (SPP, PR-PPP, or CP-PPP).

* Requires a configuration file specifying input data, processing options, and models.

### 6.2 Post-processing analysis

```bash
python algs/post_processing_gnss.py <config_performance.json>
```

* Performs additional analysis on the results (plots, error statistics, residual analysis, comparisons).

* Uses a different configuration file.

## 7. Configuration

Configuration is handled via json files. Examples for `<config_gnss.json>` and `<config_performance.json>` are provided in [configs](configs).

* [brux_2024_gnss_solver.json](configs/brux_2024_gnss_solver.json), [jfng_gnss_solver.json](configs/jfng_gnss_solver.json), [zim_gnss_solver.json](configs/zim_gnss_solver.json) for GNSS algorithms
* [brux_performance.json](configs/brux_performance.json), [jfng_performance.json](configs/jfng_performance.json), [zim_performance.json](configs/zim_performance.json) for post-processing algorithms

See the html documentation with json schemas for both configurations in [gnss_doc.html](configs/gnss_doc.html) and [performance_doc.html](configs/performance_doc.html).


## 8. Project Structure

```graphql
PyGNSSFix/
│── algs/                   # Algorithm entry points and executable scripts
│── configs/                # Configuration templates and parameter files
│── runs/                   # Generated outputs from execution runs
│── src/                    # Core source code (modules, classes, utilities)
│── docs/                   # Documentation and theory notes (algorithms, models)
│── workspace/              # Working directory for local data and experiments
│   ├── datasets/           # GNSS input datasets (RINEX, SP3, CLK, etc.)
│   ├── faults/             # Fault injection or error scenario files
│   └── geo_time_data/      # Auxiliary geophysical/temporal data (CSpice, EOP, ITRF, tropo)
│── LICENSE                 # License information
│── requirements.txt        # Python dependencies
│── README.md               # Project overview and usage instructions

```


## 9. References

The following reference books/references have been used to implement this project.

[1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck, Springer Cham, 2017

[2] ESA GNSS DATA PROCESSING, Volume I: Fundamentals and Algorithms, J. Sanz Subirana, J.M. Juan Zornoza and M. Hernández-Pajares

[3] European GNSS (GALILEO) Open Service OS, Signal-In-Space Interface Control Document (ICD). Issue 2.0, January 2021

[4] NAVSTAR GPS Space Segment/Navigation User Interfaces (IS-GPS-200). May 2021


These references are used in the documentation of some classes/functions.

## 10. License

This project is built under an [MIT License](LICENSE).






