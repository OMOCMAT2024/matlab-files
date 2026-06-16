# TOTPT Chapter 9 corrected EOM / tyre-load implementation, audited v2

This package modifies the uploaded TOTPT MATLAB/CasADi files to implement the corrected Chapter 9 roll/pitch/yaw handling model with the corresponding Chapter 9 tyre-load equations.

## Main model formulation

The state vector is

```matlab
x = [u, v, r, phi, phidot, theta, thetadot, n, xi, omega_fl, omega_fr, omega_rl, omega_rr]
```

The lifted algebraic vector is only the four tyre vertical loads:

```matlab
z = [Fz_fl, Fz_fr, Fz_rl, Fz_rr]
```

The accelerations `ax`, `ay`, `r_dot`, `phi_ddot`, and `theta_ddot` are computed explicitly from the corrected Chapter 9 EOM after tyre loads and tyre forces are known. They are not introduced as additional NLP algebraic variables.

## Corrected Chapter 9 equations implemented

With `H = h - d`, the model solves the corrected compact Chapter 9 EOM as

```text
ax = X/m
ay = Y/m
theta_ddot = (MM - H*X)/Iy
[Jx  -Jzx; -Jzx  Jz] [phi_ddot; r_dot] = [LM + H*Y; N - H*X*phi]
u_dot = ax + v*r - H*theta_ddot
v_dot = ay - u*r + H*phi_ddot
```

with

```text
LM = -kphi*phi - cphi*phidot + m*g*H*phi
MM = -ktheta*theta - ctheta*thetadot
```

The Chapter 9 tyre-load equations are enforced as algebraic equality path constraints:

```text
DeltaZ  = (-X*h + Iy*theta_ddot)/l
Z1      = m*g*a2/l + DeltaZ - Za1
Z2      = m*g*a1/l - DeltaZ - Za2
DeltaZ1 = (Y1*q1 + kphi1*phi + cphi1*phidot)/t1
DeltaZ2 = (Y2*q2 + kphi2*phi + cphi2*phidot)/t2
Fz_fl   = Z1/2 - DeltaZ1
Fz_fr   = Z1/2 + DeltaZ1
Fz_rl   = Z2/2 - DeltaZ2
Fz_rr   = Z2/2 + DeltaZ2
```

The errata-corrected vertical resultant is respected by using separate front and rear aero vertical terms `Za1` and `Za2`, not a repeated front term.

## Collocation modification

The original uploaded TOTPT structure used knot algebraic variables `Zk` but no collocation algebraic variables. This v2 package uses both:

```matlab
Zk  = SX.sym('Zk',  nz, N+1);
Zkj = SX.sym('Zkj', nz, N*d);
```

The collocation loop explicitly evaluates each collocation point one-by-one:

```matlab
[dXkj_j, Qkj_j] = f_dyn(Xkj(:,cidx), Uk(:,k+1), Zkj(:,cidx), kappa_col(cidx));
h_col_j = h_eq(Xkj(:,cidx), Uk(:,k+1), Zkj(:,cidx), kappa_col(cidx));
```

This avoids relying on implicit CasADi column-vectorization and guarantees that the collocation-state `Xkj(:,cidx)` and collocation-load `Zkj(:,cidx)` are used together in both dynamics and algebraic/path constraints.

## IPOPT-friendly changes

- The friction constraints are written without division by `Fz^2`:

```text
Fx^2/mux^2 + Fy^2/muy^2 <= Fz^2
```

- Positive lower bounds are imposed on `Fz` at knots and collocation points.
- The tyre-load loop is lifted as sparse algebraic constraints instead of being hidden inside a dense symbolic inverse.
- The model no longer relies on matrix-valued calls such as `f_dyn(Xkj(:,idx),...)`; this v2 version uses scalar/vector calls per collocation point.

## Important limitations

This code was statically audited in the ChatGPT container, but it was not executed with MATLAB/CasADi/IPOPT here because those runtime dependencies are unavailable. The roll/pitch numerical parameters in `veh_params.m` are placeholders and must be replaced by identified vehicle data before using the model for engineering conclusions.
