clear, clc
% addpath('D:\App\adi-3.7.2-windows64-matlab2018b')
addpath('D:\casadi_3_7_2')

veh = my_params();
track_name = char("Catalunya");

% Rosberg DRS/KERS data used to build the distance-varying schedules.
% Change this path only if the CSV is stored in another folder.
drs_data_file = fullfile(pwd, ...
    'rosberg_2012_china_distance_speed_drs_kers_1m(4).csv');

% Boundary-condition mode for the lap optimization.
% false: retain the existing non-periodic start/end boundary conditions.
% true:  impose Xk(:,end) = Xk(:,1) for all vehicle states.
use_periodic_state_bc = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath('D:\App\my_oc_test\my_shanghai_track')
addpath('D:\Applications\TOTPT-main\oc\my_tracks\Shanghai')
load C_true_project
load nl
load nr
load OmegaB
load Theta
load S_all_points
load P_L
load P_R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

curvature = OmegaB(:, 3);
% Theta = wrapToPi(Theta);
% curvature = OmegaB(:, 3) * 0;

track_matrix = [C_true_project(:, 1), C_true_project(:, 2), (-1)*nr(:, 1), nl(:, 1), Theta, curvature];

P_L(end, :) = [];
P_R(end, :) = [];

S_all_points(end, :) = [];

track_matrix(end, :) = []; % do this because the raw data's first and last points overlap

figure();plot(track_matrix(:,1), track_matrix(:,2), '-o');hold on;grid on;plot(track_matrix(1,1), track_matrix(1,2), '-s');plot(track_matrix(end,1), track_matrix(end,2), '-*');axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extend the starting line backwards a bit

% num_extend_points = 2; % 2
%
% track_matrix_points_extended = [track_matrix((length(track_matrix(:,1))-num_extend_points):end, :); track_matrix];

% S_all_points_extended = [S_all_points((length(S_all_points(:,1))-num_extend_points):end, :); S_all_points];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N_query_points = 1884; % this must be an integer (4000) 945 1884 1884*3
% N_query_points = 2046/2; % this was used for circle drifting
% N_query_points = 3930/2; % this was used for figure-eight drifting
N_query_points = 5444/2; % this was used for Shanghai
[track_interp_here, s_interp_here] = calc_track_interp(track_matrix, S_all_points, N_query_points);



[P_L_interp_here, ~] = calc_track_interp(P_L, S_all_points, N_query_points);
[P_R_interp_here, ~] = calc_track_interp(P_R, S_all_points, N_query_points);



% Build piecewise-constant Cd and KERS schedules from Rosberg's measured
% DRS/KERS states. The source distance is scaled to the optimization track
% length so the transitions remain aligned when the two lap-length
% representations differ slightly.
if ~isfile(drs_data_file)
    error('my_oc:MissingDRSFile', ...
        'DRS/KERS data file not found: %s', drs_data_file);
end

drs_data = readtable(drs_data_file);
required_drs_variables = ...
    {'distance_official_length_scaled_m','drs_on','kers_active'};
if ~all(ismember(required_drs_variables, drs_data.Properties.VariableNames))
    error('my_oc:InvalidDRSFile', ...
        ['DRS/KERS CSV must contain distance_official_length_scaled_m, ' ...
         'drs_on and kers_active.']);
end

drs_distance_m = drs_data.distance_official_length_scaled_m;
drs_on_source = drs_data.drs_on;
kers_on_source = drs_data.kers_active;

if any(~isfinite(drs_distance_m)) || any(diff(drs_distance_m) <= 0)
    error('my_oc:InvalidDRSDistance', ...
        'DRS/KERS distance data must be finite and strictly increasing.');
end
if any(~ismember(drs_on_source, [0 1]))
    error('my_oc:InvalidDRSState', 'DRS state data must contain only 0 or 1.');
end
if any(~ismember(kers_on_source, [0 1]))
    error('my_oc:InvalidKERSState', ...
        'KERS state data must contain only 0 or 1.');
end
if s_interp_here(end) <= 0 || drs_distance_m(end) <= 0
    error('my_oc:InvalidLapLength', ...
        'Track and DRS/KERS lap lengths must be positive.');
end

% Identify the final contiguous KERS-on interval directly from the binary
% source data. This avoids hard-coding its distance limits.
kers_start_idx = find(diff([0; double(kers_on_source(:))]) == 1);
kers_end_idx = find(diff([double(kers_on_source(:)); 0]) == -1);
if isempty(kers_start_idx) || numel(kers_start_idx) ~= numel(kers_end_idx)
    error('my_oc:InvalidKERSIntervals', ...
        'Could not identify valid contiguous KERS-on intervals.');
end
kers_last_on_source = zeros(size(kers_on_source));
kers_last_on_source(kers_start_idx(end):kers_end_idx(end)) = 1;

drs_query_distance_m = ...
    s_interp_here ./ s_interp_here(end) .* drs_distance_m(end);
drs_on_track = interp1(drs_distance_m, double(drs_on_source), ...
    drs_query_distance_m, 'previous', 'extrap');
kers_on_track = interp1(drs_distance_m, double(kers_on_source), ...
    drs_query_distance_m, 'previous', 'extrap');
kers_last_on_track = interp1(drs_distance_m, ...
    double(kers_last_on_source), drs_query_distance_m, ...
    'previous', 'extrap');

% Restore exact binary KERS schedules after interpolation.
kers_on_track = double(kers_on_track >= 0.5);
kers_last_on_track = double(kers_last_on_track >= 0.5);
if any(kers_last_on_track > kers_on_track)
    error('my_oc:InvalidLastKERSInterval', ...
        'The final KERS interval must be a subset of the KERS-on schedule.');
end

Cd_track = veh.Cd_nominal * ones(size(s_interp_here));
Cd_track(drs_on_track >= 0.5) = veh.Cd_drs;

track_interp_table = array2table([track_interp_here(:, 1:4) ...
    track_interp_here(:, 5) track_interp_here(:, 6) s_interp_here ...
    Cd_track kers_on_track kers_last_on_track]);
track_interp_table.Properties.VariableNames = ...
    {'x','y','wr','wl','Theta','curvature','s','Cd', ...
     'kers_on','kers_last_on'};

if ~exist(fullfile(pwd,'my_smooth'), 'dir')
mkdir('my_smooth')
end
writetable(track_interp_table,['my_smooth/' track_name '_smooth' '.csv'])







%%
run('my_nlp.m')

my_plot

save('result', 'sol', 'track')
