# Simple Time Simulation v3 checked

This package contains one MATLAB file:

- `run_oc_simple_time_sim_arbitrary.m`

Changes relative to v3:

1. Keeps the true `ode15i` DAE formulation for `[x; z]`.
2. Leaves fixed-step RK4 state bounds unclamped by default so RK4, ode15s, and ode15i solve the same simulation model. You can turn state clamping back on with `opts.enforce_state_bounds_fixed = true`.
3. Adds a diagnostic `sim_ode15i_dae.decic_residual_inf` to check the consistency of the DAE initial condition.

Recommended test:

```matlab
clear functions
rehash
opts = struct();
opts.t_end = 3;
opts.fixed_dt = 1e-3;
opts.do_plot = true;
[sim_ode15s, sim_fixed, sim_ode15i_dae, model, ctrl] = run_oc_simple_time_sim_arbitrary('mild_sine', opts);
max(abs(sim_ode15s.h_load(:)))
max(abs(sim_fixed.h_load(:)))
max(abs(sim_ode15i_dae.h_load(:)))
sim_ode15i_dae.decic_residual_inf
```
