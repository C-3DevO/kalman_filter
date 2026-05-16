# Kalman Filter for Time-Varying Channel Estimation

This project implements a Kalman filter for estimating the coefficients of a time-varying FIR communication channel using MATLAB.

The simulation investigates recursive channel estimation under noisy conditions and studies how the Kalman filter tracks slowly and rapidly varying channel dynamics.

---

## Overview

The communication channel is modeled as a time-varying FIR filter whose tap coefficients evolve according to a Gauss–Markov process.

The observation model is given by:

```math
x[n] = v^T[n]h[n] + w[n]
```

where:

- `v[n]` is the known transmitted input signal
- `h[n]` is the unknown channel coefficient vector
- `w[n]` is additive white Gaussian noise (AWGN)

The channel evolution follows the state-space model:

```math
h[n] = Ah[n-1] + u[n]
```

The Kalman filter recursively estimates the channel coefficients using prediction and correction steps.

---

## Contents

- MATLAB implementation of Kalman filtering
- Time-varying FIR channel modeling
- Gauss–Markov state-space representation
- Recursive MMSE channel estimation
- Kalman gain analysis
- Estimation error covariance analysis
- Tracking performance under different channel dynamics

---

## Simulations Performed

### Task 1 — Slow-Fading Channel Estimation

Simulation using:

- Diagonal state transition matrix
- Small process noise covariance
- Slowly varying channel taps

The Kalman filter successfully converges and tracks the channel coefficients with low estimation error.

### Task 2 — Increased Process Noise

The process noise covariance `Q` is increased to create faster channel variations.

Observed effects:

- Faster channel fluctuations
- Larger Kalman gains
- Higher steady-state estimation error
- Reduced tracking smoothness

### Task 3 — Coupled Channel Dynamics

Off-diagonal elements are introduced in the state transition matrix:

```math
A =
\begin{bmatrix}
0.99 & 0.01 \\
0.01 & 0.999
\end{bmatrix}
```

This introduces correlation between channel taps and produces more complex channel dynamics.

---

## Kalman Filter Equations

### Prediction Step

```math
\hat{h}[n|n-1] = A\hat{h}[n-1|n-1]
```

```math
M[n|n-1] = AM[n-1|n-1]A^T + Q
```

### Kalman Gain

```math
K[n] =
\frac{M[n|n-1]v[n]}
{\sigma^2 + v^T[n]M[n|n-1]v[n]}
```

### Update Step

```math
\hat{h}[n|n] =
\hat{h}[n|n-1]
+ K[n]\left(x[n]-v^T[n]\hat{h}[n|n-1]\right)
```

---

## Observations

- The Kalman filter converges rapidly after initialization
- Estimated channel taps closely follow the true channel coefficients
- Larger process noise increases tracking difficulty
- Coupled channel dynamics lead to more complex estimation behavior
- Kalman gains decrease as uncertainty reduces
- Estimation MSE stabilizes after convergence

---

## References

1. S. M. Kay, *Fundamentals of Statistical Signal Processing: Estimation Theory*, Prentice Hall, 1993.
