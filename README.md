# Recursive State Estimation for Dynamic Systems

## 📍 Overview

This repository contains a complete MATLAB implementation for comparing different state estimation methods used in tracking the 2D position and velocity of an autonomous vehicle on a racetrack. The methods include various Kalman filters and particle filtering techniques.

The dataset and task were provided as part of the course *Recursive State Estimation for Dynamic Systems* at the Geodätisches Institut, Leibniz Universität Hannover (SoSe 2024).

## 🚗 Scenario
![image](https://github.com/user-attachments/assets/be103192-2a45-464f-bdc1-49b7eb97776c)

An autonomous test vehicle moves along a predefined racetrack, tracked by two total stations (TS1 & TS2). Each station provides:
- Distance `d_TSi(k)` and
- Horizontal angle `α_TSi(k)` at each time step `k`.

A section of the track has a known constant curve radius `R = 9.00 m`. This geometrical constraint is used in some of the filtering approaches.

## 📦 State and Observation Models

### 🔧 State vector `x(k)`:
```
x(k) = [x(k), y(k), v_x(k), v_y(k)]ᵀ
```
2D position and velocity in meters and m/s.

### 📡 Observation vector `l(k)`:
```
l(k) = [d_TS1(k), α_TS1(k), d_TS2(k), α_TS2(k)]
```

Known parameters:
- TS1 position: [135.54 m, 98.79 m]
- TS2 position: [110.00 m, 90.00 m]

## 📐 Implemented Filters

- [x] Extended Kalman Filter (EKF)
- [x] EKF with known curve constraint
- [x] Iterated EKF (IEKF)
- [x] Unscented Kalman Filter (UKF)
- [x] Ensemble Kalman Filter (EnKF)
- [x] Particle Filter (PF)

Each filter is modularly implemented in separate MATLAB functions.

## 📊 Tasks

- Build system and measurement models (`ffun.m` and `hfun.m`)
- Implement multiple filters (`filter_*.m`)
- Perform RMSE-based comparison over time
- Conduct Monte Carlo simulations for parameter sensitivity:
  - Iterations (IEKF)
  - Samples (EnKF)
  - Particles (PF)
- Analyze effect of sampling rate and noise parameters

## 📁 File Structure

```
main.m                   % Main entry for data loading and execution
ffun.m                   % System model
hfun.m                   % Measurement model
filter_ekf.m             % EKF implementation
filter_iekf.m            % IEKF implementation
filter_ukf.m             % UKF implementation
filter_enkf.m            % EnKF implementation
filter_pf.m              % Particle filter implementation
...
```


## 🛠 Notes

- Use `atan2()` for angle computation to ensure correct quadrant detection.
- Pay attention to unit conversions (e.g., mm vs m, gon vs rad).
- Compare filter output with ground truth for plausibility.
- Use symbolic differentiation (e.g. WolframAlpha) for Jacobians in EKF.

## 📬 Contact

M.Sc. Tingde Liu 
E-Mail: tingde.liu.luh@gmail.com 


---

> This repository is part of coursework under Prof. Dr.-Ing. Ingo Neumann, Leibniz Universität Hannover.

