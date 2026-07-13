DISTANCE-BASED DRS + KERS IMPLEMENTATION — RE-AUDITED
=====================================================

Purpose
-------
This package retains the distance-based DRS implementation and adds
Rosberg's distance-based KERS ON/OFF schedule to the existing direct optimal
control model. The four core MATLAB files are unchanged from the preceding
DRS+KERS package because the re-audit found no core-code correction needed.
The validation and audit documents have been strengthened and corrected.

Distance-based DRS
------------------
DRS off:
    Cd = veh.Cd_nominal = 0.9

DRS on:
    Cd = veh.Cd_drs = 0.8

The drag force in my_model.m is:
    f_drag = 0.5 * veh.rho * Cd * veh.A * ux.^2

Distance-based KERS
-------------------
The Rosberg CSV contains four KERS-on intervals. my_oc.m identifies the final
contiguous interval automatically, without hard-coded distance limits.

my_model.m uses:
    pv = [kappa; Cd; kers_on; kers_last_on]

The available rear-drivetrain power ceiling in my_nlp.m is:

KERS off:
    veh.power_motor_max

Earlier KERS intervals:
    veh.power_motor_max + veh.kers_power_add

Final KERS interval:
    veh.power_motor_max + veh.kers_power_add_last

Tuning the final KERS interval
------------------------------
Edit only this line in my_params.m:

    veh.kers_power_add_last = 60 * 10^3;

The 2012 FIA maximum KERS power was 60 kW. A value above 60 kW is therefore
an empirical model-calibration setting rather than a regulation-faithful
KERS power setting.

400 kJ/lap distinction
----------------------
The 2012 FIA regulations separately limited energy released from KERS to
400 kJ per lap. The current OCP implements the requested distance-based
increase in the available-power ceiling. It does not introduce a KERS-energy
state or an integral energy constraint.

The extracted KERS-on durations total approximately 6.9805 s. If the added
power ceiling were saturated continuously at 60 kW throughout all four
windows, the maximum possible added energy would be approximately 418.83 kJ.
This does not mean the solved trajectory necessarily uses 418.83 kJ: the
optimizer is permitted, not forced, to saturate the added power ceiling.
The validator reports this distinction explicitly.

Files changed relative to the user's originals
----------------------------------------------
my_params.m:
    Adds only DRS aerodynamic parameters and the two KERS power parameters.

my_oc.m:
    Reads the combined DRS/KERS CSV, maps the binary schedules to track
    distance, and adds Cd, kers_on and kers_last_on to the track table.

my_model.m:
    Expands pv and uses distance-varying Cd in drag.

my_nlp.m:
    Passes the four-row pv everywhere and increases the existing power
    ceiling only in the KERS-on windows.

my_plot.m:
    Byte-for-byte unchanged.

No state, control, algebraic variable, tyre model, suspension model, lift or
downforce model, torque distribution, bound, scaling, initial guess,
objective, collocation degree, periodic boundary condition, IPOPT option, or
plotting logic was changed.

Validation
----------
Place all package files in one working folder and run:

    validate_drs_kers_setup

Then run:

    my_oc

The validator intentionally issues a warning when the maximum possible
full-saturation added energy exceeds 400 kJ. This is a modelling-scope
warning, not a MATLAB-code failure.
