# Simple time simulation package v2

This package contains a single MATLAB simulator:

- `run_oc_simple_time_sim_arbitrary.m`

It builds CasADi functions directly from `my_model.m` and `my_params.m`, generates arbitrary drive torque and steering controls, sets brake torque to zero, applies the rear-drive power limit online, and runs both:

1. `ode15s` variable-step stiff simulation of the reduced ODE `xdot = f(x,u,z*)`, where `z*` solves `h_load = 0`.
2. fixed-step RK4 with `h_load = 0` solved at every RK stage.

Version v2 uses the analytic CasADi Jacobian `dh_load/dz` for the four tyre-load algebraic equations. This replaces the finite-difference Jacobian in v1 and should make RK4 much faster.

Example:

```matlab
clear functions
rehash
addpath('D:\casadi_3_7_2')
addpath('D:\Applications\TOTPT-main\oc')

opts = struct();
opts.t_end = 3;
opts.fixed_dt = 1e-3;
opts.do_plot = true;

[sim_ode15s, sim_fixed, model, ctrl] = run_oc_simple_time_sim_arbitrary('mild_sine', opts);
```

Check:

```matlab
max(abs(sim_ode15s.h_load(:)))
max(abs(sim_fixed.h_load(:)))
```
