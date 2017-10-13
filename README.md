# CarND-Controls-MPC
---

The purpose of this project is to develop a model predictive controller (MPC) with additional latency 100ms between actuator commands to steer a car. The vehicle model used is a kinematic model. It neglects many dynamical effects such as friction and tire sliding. 

## Kinematic Vehicle Model
The model uses following equations:
```
  // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
  // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
  // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
  // v_[t+1] = v[t] + a[t] * dt
  // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
  // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt

```
* `x,y`: position of the car
* `psi`: heading direction
* `v`: velocity 
* `cte`: cross-track error 
* `epsi`: orientation error 
* `Lf`: distance between the center of mass of the vehicle and the front wheels 

## Timestep Length and Elapsed Duration (N & dt)

The time frame `T=N*dt` defines the update rate. Higher rate (lower `T`) lead to more responsive controlers, but can suffer from instabilities. Lower rate generally leads to smoother controls. For a given rate, shorter time steps `dt` imply finer controls but also require a larger MPC problem to be solved. 

In the projcet, `N=9` and `dt=0.12` are set.

## Polynomial Fitting and MPC Preprocessing

The waypoints are preprocessed, they are transformed to the vehicle's perspective. This simplifies the process to fit a polynomial to the waypoints. Threrefore, the vehicle's x and y coordinates are now at the origin (0, 0) and the orientation angle is also 0.

## Model Predictive Control with Latency

The original kinematic equations depend on the actuations from the previous timestep, but with a delay of 100ms (which happens to be the timestep interval) the actuations are applied another timestep later, so the equations have been altered to account for this (MPC.cpp lines 95-97). 
The cost functions, such as punishing CTE, epsi, difference between velocity and a reference velocity, delta, acceleration, change in delta, and change in acceleration are fine tuned based on experiment results. 

---
# Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.
