% left/right boundaries as Nx2 arrays:
% boundL = [xL yL]; boundR = [xR yR];

opts = struct();
opts.stepsize_prep = 1.0;
opts.stepsize_opt  = 3.0;   % like TUMFTM default :contentReference[oaicite:7]{index=7}
opts.stepsize_out  = 2.0;
opts.width_opt     = 3.4;   % vehicle width + safety (TUMFTM default) :contentReference[oaicite:8]{index=8}
opts.curvlim       = 0.12;  % steering curvature limit (paper/repo default) :contentReference[oaicite:9]{index=9}
opts.use_curv_constr = true;

out = mincurv_raceline_from_bounds(boundL, boundR, opts);

figure; hold on; axis equal; grid on;
plot(boundL(:,1), boundL(:,2), 'k-');
plot(boundR(:,1), boundR(:,2), 'k-');
plot(out.raceline_xy(:,1), out.raceline_xy(:,2), 'r-', 'LineWidth', 2);
legend('Left bound','Right bound','Min-curv raceline');