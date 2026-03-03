% clc; clear; close all;
% 
% fname = "edge_loop_00.csv";   % put in current folder, or use full path
% T = readtable(fname);         % reads header x,y,z automatically
% 
% x = T{:,1};
% y = T{:,2};
% z = T{:,3};
% 
% figure;
% plot3(x, y, z, '-', 'LineWidth', 1.5); hold on;
% scatter3(x(1), y(1), z(1), 60, 'filled');           % start point
% scatter3(x(end), y(end), z(end), 60, 'filled');     % end point
% grid on; axis equal; view(3);
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title('3D trajectory from CSV');





clc; clear; close all;

% Folder that contains the CSVs (use pwd for current folder)
folder = pwd;   % or: folder = "C:\path\to\your\csvs";

files = dir(fullfile(folder, "edge_loop_*.csv"));
assert(~isempty(files), "No files matched edge_loop_*.csv in: %s", folder);

figure; hold on; grid on; axis equal; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('All edge\_loop\_*.csv');

for k = 1:numel(files)
    fname = fullfile(files(k).folder, files(k).name);

    T = readtable(fname);     % expects header x,y,z (works even if header differs)
    x = T{:,1}; y = T{:,2}; z = T{:,3};

    plot3(x, y, z, '-', 'LineWidth', 1.2);

    % Optional: mark start/end points
    scatter3(x(1),   y(1),   z(1), 25, 'filled');
    scatter3(x(end), y(end), z(end), 25, 'filled');
end

legend({files.name}, 'Interpreter','none', 'Location','bestoutside'); % optional (can be long)