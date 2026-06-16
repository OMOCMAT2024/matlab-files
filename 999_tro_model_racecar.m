import casadi.*

veh = veh_params();

%% Vehicle Model - state variables
% Corrected Chapter 9 model with velocity at the vehicle invariant point M.
% States are scaled.  Physical variables are recovered by multiplying by x_s.

nx = 13;

% VIP longitudinal velocity u [m/s]
ux_n = SX.sym('ux_n');
ux_s = veh.ux_s;
ux = ux_s * ux_n;

% VIP lateral velocity v [m/s]
uy_n = SX.sym('uy_n');
uy_s = veh.uy_s;
uy = uy_s * uy_n;

% yaw rate r [rad/s]
r_n = SX.sym('r_n');
r_s = veh.r_s;
r = r_s * r_n;

% roll angle phi [rad]
phi_n = SX.sym('phi_n');
phi_s = veh.phi_s;
phi = phi_s * phi_n;

% roll rate phidot [rad/s]
phidot_n = SX.sym('phidot_n');
phidot_s = veh.phidot_s;
phidot = phidot_s * phidot_n;

% pitch angle theta [rad]
theta_n = SX.sym('theta_n');
theta_s = veh.theta_s;
theta = theta_s * theta_n;

% pitch rate thetadot [rad/s]
thetadot_n = SX.sym('thetadot_n');
thetadot_s = veh.thetadot_s;
thetadot = thetadot_s * thetadot_n;

% lateral distance to centreline [m] - left of centreline => n > 0
n_n = SX.sym('n_n');
n_s = veh.n_s;
n = n_s * n_n;

% course angle to centreline tangent direction [rad]
xi_n = SX.sym('xi_n');
xi_s = veh.xi_s;
xi = xi_s * xi_n;

omega_s = veh.omega_s;

omega_fl_n = SX.sym('omega_fl_n');
omega_fr_n = SX.sym('omega_fr_n');
omega_rl_n = SX.sym('omega_rl_n');
omega_rr_n = SX.sym('omega_rr_n');
omega_fl = omega_s * omega_fl_n;
omega_fr = omega_s * omega_fr_n;
omega_rl = omega_s * omega_rl_n;
omega_rr = omega_s * omega_rr_n;

x = [ux_n; uy_n; r_n; phi_n; phidot_n; theta_n; thetadot_n; n_n; xi_n; ...
     omega_fl_n; omega_fr_n; omega_rl_n; omega_rr_n];

x_s = [ux_s; uy_s; r_s; phi_s; phidot_s; theta_s; thetadot_s; n_s; xi_s; ...
       omega_s; omega_s; omega_s; omega_s];

x_min = [veh.ux_min; veh.uy_min; veh.r_min; veh.phi_min; veh.phidot_min; ...
         veh.theta_min; veh.thetadot_min; veh.n_min; veh.xi_min; ...
         veh.omega_min; veh.omega_min; veh.omega_min; veh.omega_min] ./ x_s;
x_max = [veh.ux_max; veh.uy_max; veh.r_max; veh.phi_max; veh.phidot_max; ...
         veh.theta_max; veh.thetadot_max; veh.n_max; veh.xi_max; ...
         veh.omega_max; veh.omega_max; veh.omega_max; veh.omega_max] ./ x_s;

%% Vehicle model - control variables
nu = 3;

Tt_n = SX.sym('Tt_n');
Tt_s = veh.Tt_s;
Tt = Tt_s * Tt_n;

% Braking torque is negative as in the original TOTPT convention.
Tb_n = SX.sym('Tb_n');
Tb_s = veh.Tb_s;
Tb = Tb_s * Tb_n;

delta_n = SX.sym('delta_n');
delta_s = veh.delta_s;
delta = delta_s * delta_n;

u = [Tt_n; Tb_n; delta_n];
u_s = [Tt_s; Tb_s; delta_s];

u_min = [veh.Tt_min; veh.Tb_min; veh.delta_min] ./ u_s;
u_max = [veh.Tt_max; veh.Tb_max; veh.delta_max] ./ u_s;

duk_lb = [veh.Tt_dot_min; veh.Tb_dot_min; veh.delta_dot_min] ./ u_s;
duk_ub = [veh.Tt_dot_max; veh.Tb_dot_max; veh.delta_dot_max] ./ u_s;

% Regularisation factors for knot-to-knot rates used in tro_code.m.
rdu = [1; 1; 10];
rdy = zeros(4,1);

%% Algebraic / auxiliary variables
% Lift only the four tyre vertical loads.  This breaks the tyre-load algebraic
% loop without lifting ax, ay, r_dot, phi_ddot or theta_ddot.
nz = 4;

fzfl_n = SX.sym('fzfl_n');
fzfr_n = SX.sym('fzfr_n');
fzrl_n = SX.sym('fzrl_n');
fzrr_n = SX.sym('fzrr_n');

Fz_s = veh.Fz_s;
fzfl = Fz_s * fzfl_n;
fzfr = Fz_s * fzfr_n;
fzrl = Fz_s * fzrl_n;
fzrr = Fz_s * fzrr_n;

z = [fzfl_n; fzfr_n; fzrl_n; fzrr_n];
z_s = [Fz_s; Fz_s; Fz_s; Fz_s];
z_min = [veh.Fz_min; veh.Fz_min; veh.Fz_min; veh.Fz_min] ./ z_s;
z_max = [veh.Fz_max; veh.Fz_max; veh.Fz_max; veh.Fz_max] ./ z_s;

%% Variable parameter
% kappa > 0 for left turns.
kappa = SX.sym('kappa');
pv = kappa;

%% Vehicle model equations
cosd = cos(delta);
sind = sin(delta);

% Aerodynamics.  With this sign convention, f_lift < 0 is downforce.
f_drag = 0.5 * veh.drag_coeff * ux.^2;
f_lift = 0.5 * veh.lift_coeff * ux.^2;
Za1 = veh.aero_front_frac * f_lift;
Za2 = (1 - veh.aero_front_frac) * f_lift;

% Tyre contact-point velocities in wheel frames.  Roll and pitch affect the
% Chapter 9 body EOM and vertical loads; they do not directly modify the tyre
% constitutive equations in this small-angle rigid-tyre model.
vxfl_raw = (ux - 0.5*veh.t1*r)*cosd + (uy + veh.lf*r)*sind;
vxfr_raw = (ux + 0.5*veh.t1*r)*cosd + (uy + veh.lf*r)*sind;
vxrl_raw =  ux - 0.5*veh.t2*r;
vxrr_raw =  ux + 0.5*veh.t2*r;

vyfl = (uy + veh.lf*r)*cosd - (ux - 0.5*veh.t1*r)*sind;
vyfr = (uy + veh.lf*r)*cosd - (ux + 0.5*veh.t1*r)*sind;
vyrl =  uy - veh.lr*r;
vyrr =  uy - veh.lr*r;

% Smoothly protect the tyre-slip divisions from zero wheel-frame forward speed.
% In normal race operation vx*_raw remains positive and this equals vx*_raw to
% numerical precision.  The small guard is only for IPOPT robustness.
vx_eps = veh.vx_eps;
vxfl = regularized_positive_den(vxfl_raw, vx_eps);
vxfr = regularized_positive_den(vxfr_raw, vx_eps);
vxrl = regularized_positive_den(vxrl_raw, vx_eps);
vxrr = regularized_positive_den(vxrr_raw, vx_eps);

SXfl = (veh.rw*omega_fl - vxfl_raw) ./ vxfl;
SXfr = (veh.rw*omega_fr - vxfr_raw) ./ vxfr;
SXrl = (veh.rw*omega_rl - vxrl_raw) ./ vxrl;
SXrr = (veh.rw*omega_rr - vxrr_raw) ./ vxrr;

syfl = vyfl ./ vxfl;
syfr = vyfr ./ vxfr;
syrl = vyrl ./ vxrl;
syrr = vyrr ./ vxrr;

sfl = sqrt(SXfl.^2 + syfl.^2);
sfr = sqrt(SXfr.^2 + syfr.^2);
srl = sqrt(SXrl.^2 + syrl.^2);
srr = sqrt(SXrr.^2 + syrr.^2);

slip_eps = veh.slip_eps;
sfl_d = sqrt(sfl.^2 + slip_eps.^2);
sfr_d = sqrt(sfr.^2 + slip_eps.^2);
srl_d = sqrt(srl.^2 + slip_eps.^2);
srr_d = sqrt(srr.^2 + slip_eps.^2);

fxfl = (SXfl ./ sfl_d) .* MF_Fx(veh, sfl, fzfl);
fxfr = (SXfr ./ sfr_d) .* MF_Fx(veh, sfr, fzfr);
fxrl = (SXrl ./ srl_d) .* MF_Fx(veh, srl, fzrl);
fxrr = (SXrr ./ srr_d) .* MF_Fx(veh, srr, fzrr);
fxf = fxfl + fxfr;
fxr = fxrl + fxrr;

fyfl = (syfl ./ sfl_d) .* MF_Fy(veh, sfl, fzfl);
fyfr = (syfr ./ sfr_d) .* MF_Fy(veh, sfr, fzfr);
fyrl = (syrl ./ srl_d) .* MF_Fy(veh, srl, fzrl);
fyrr = (syrr ./ srr_d) .* MF_Fy(veh, srr, fzrr);
fyf = fyfl + fyfr;
fyr = fyrl + fyrr;

mufl = sqrt(fxfl.^2 + fyfl.^2) ./ fzfl;
mufr = sqrt(fxfr.^2 + fyfr.^2) ./ fzfr;
murl = sqrt(fxrl.^2 + fyrl.^2) ./ fzrl;
murr = sqrt(fxrr.^2 + fyrr.^2) ./ fzrr;

%% Wheel torque allocation
fzf = fzfl + fzfr;
fzr = fzrl + fzrr;
dfzf = 0.5 * (fzfr - fzfl);
dfzr = 0.5 * (fzrr - fzrl);

Tf = veh.kb * Tb;
Tfl = 0.5 * Tf * (1 - dfzf ./ fzf);
Tfr = 0.5 * Tf * (1 + dfzf ./ fzf);

Tr = Tt + (1 - veh.kb) * Tb;
Trl = 0.5 * Tr * (1 - dfzr ./ fzr);
Trr = 0.5 * Tr * (1 + dfzr ./ fzr);

%% Corrected Chapter 9 force and moment resultants
X1 = fxf*cosd - fyf*sind;
X2 = fxr;
Y1 = fyf*cosd + fxf*sind;
Y2 = fyr;

X = X1 + X2 - f_drag;
Y = Y1 + Y2;

N = veh.lf*Y1 - veh.lr*Y2 ...
    + 0.5*veh.t1*((fxfr - fxfl)*cosd + (fyfl - fyfr)*sind) ...
    + 0.5*veh.t2*(fxrr - fxrl);

H = veh.H;
LM = -veh.kphi*phi - veh.cphi*phidot + veh.m*veh.g*H*phi;
MM = -veh.ktheta*theta - veh.ctheta*thetadot;

% Corrected Chapter 9 compact EOM, solved explicitly after X,Y,N are known.
ax = X ./ veh.m;
ay = Y ./ veh.m;

theta_ddot = (MM - H*X) ./ veh.Iy;

Djr = veh.Ix*veh.Iz - veh.Jzx^2;
Aphi = LM + H*Y;
Ar = N - H*X*phi;
phi_ddot = (veh.Iz*Aphi + veh.Jzx*Ar) ./ Djr;
r_dot = (veh.Jzx*Aphi + veh.Ix*Ar) ./ Djr;

ux_dot = ax + uy*r - H*theta_ddot;
uy_dot = ay - ux*r + H*phi_ddot;

%% Corrected Chapter 9 tyre-load equations, enforced as path constraints
Z10 = veh.m*veh.g*veh.lr/veh.l;
Z20 = veh.m*veh.g*veh.lf/veh.l;

DeltaZ = (-X*veh.h + veh.Iy*theta_ddot) ./ veh.l;
Z1 = Z10 + DeltaZ - Za1;
Z2 = Z20 - DeltaZ - Za2;

DeltaZ1 = (Y1*veh.q1 + veh.kphi1*phi + veh.cphi1*phidot) ./ veh.t1;
DeltaZ2 = (Y2*veh.q2 + veh.kphi2*phi + veh.cphi2*phidot) ./ veh.t2;

fzfl_ch9 = 0.5*Z1 - DeltaZ1;
fzfr_ch9 = 0.5*Z1 + DeltaZ1;
fzrl_ch9 = 0.5*Z2 - DeltaZ2;
fzrr_ch9 = 0.5*Z2 + DeltaZ2;

load_res = [fzfl - fzfl_ch9;
            fzfr - fzfr_ch9;
            fzrl - fzrl_ch9;
            fzrr - fzrr_ch9];

h_load = load_res ./ Fz_s;

%% Change of independent variable and state derivatives
% sf = dt/ds.
sf_den_raw = ux*cos(xi) - uy*sin(xi);
sf_den = regularized_positive_den(sf_den_raw, vx_eps);
sf = (1 - n*kappa) ./ sf_den;

n_dot = ux*sin(xi) + uy*cos(xi);
xi_dot = r - kappa ./ sf;

xdot = [
    ux_dot;
    uy_dot;
    r_dot;
    phidot;
    phi_ddot;
    thetadot;
    theta_ddot;
    n_dot;
    xi_dot;
    (Tfl - fxfl*veh.rw) ./ veh.Iw;
    (Tfr - fxfr*veh.rw) ./ veh.Iw;
    (Trl - fxrl*veh.rw) ./ veh.Iw;
    (Trr - fxrr*veh.rw) ./ veh.Iw;
    ] ./ x_s;

dx = sf .* xdot;

%% Output variables
y_aero = [f_drag; f_lift; Za1; Za2];
y_tyre = [fxfl; fxfr; fxrl; fxrr; fyfl; fyfr; fyrl; fyrr; fzfl; fzfr; fzrl; fzrr];
y_slip = [SXfl; SXfr; SXrl; SXrr; atan(syfl); atan(syfr); atan(syrl); atan(syrr)];
y_acc = [ax; ay; ux_dot; uy_dot; r_dot; phi_ddot; theta_ddot];
y_torque = [Tfl; Tfr; Trl; Trr];
y_mu = [mufl; mufr; murl; murr];
y_pwr = [Tfl*omega_fl; Tfr*omega_fr; Trl*omega_rl; Trr*omega_rr];
y_xdot = xdot;
y_load = [h_load; fzfl_ch9; fzfr_ch9; fzrl_ch9; fzrr_ch9; DeltaZ; DeltaZ1; DeltaZ2];

function y = regularized_positive_den(x, eps_val)
% Smooth guard for denominators that are expected to stay positive in the OCP.
y = x + eps_val;
end

function Fx = MF_Fx(veh,lambda,Fz)
Fx = (veh.d1x.*Fz + veh.d2x) .* sin(veh.Cx .* atan(veh.Bx .* lambda));
end

function Fy = MF_Fy(veh,alpha,Fz)
Fy = -(veh.d1y.*Fz + veh.d2y) .* sin(veh.Cy .* atan(veh.By .* alpha));
end
