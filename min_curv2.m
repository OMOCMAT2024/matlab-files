function out = mincurv_raceline_from_bounds(boundL_xy, boundR_xy, opts)
% Minimum-curvature raceline (QP + optional IQP loop), consistent with:
% - Heilmeier et al. (VSD 2019): r = p + alpha*n, alpha bounds, curvature QP, IQP loop
% - TUMFTM repo convention: normal vectors point to the RIGHT in driving direction
%
% INPUT
%   boundL_xy: [N x 2] left boundary polyline (closed or open)
%   boundR_xy: [N x 2] right boundary polyline (closed or open)
%   opts: struct (optional)
%
% OUTPUT
%   out.reftrack      [N x 4] [x y w_tr_right w_tr_left]
%   out.raceline_xy   [N x 2]
%   out.alpha         [N x 1]
%   out.kappa         [N x 1]
%   out.s             [N x 1]
%   out.debug         struct

    if nargin < 3, opts = struct(); end
    opts = set_default_opts(opts);

    % --- 1) Ensure closed ---
    boundL_xy = ensure_closed(boundL_xy);
    boundR_xy = ensure_closed(boundR_xy);

    % --- 2) Resample both bounds to a common station grid (simple, assumes roughly corresponding loops) ---
    [bL_prep, sL, LL] = resample_closed_polyline(boundL_xy, opts.stepsize_prep);
    [bR_prep, sR, LR] = resample_closed_polyline(boundR_xy, opts.stepsize_prep);

    L = min(LL, LR);
    s_common = (0:opts.stepsize_prep:(L - opts.stepsize_prep))';
    bL = interp_closed(bL_prep, sL, s_common);
    bR = interp_closed(bR_prep, sR, s_common);

    % --- 3) Initial reference line: midline ---
    pref = 0.5 * (bL + bR);
    if opts.do_smooth
        pref = smooth_closed(pref, opts.smooth_window);
    end

    % --- 4) Resample reference to optimization step (typically 3m) ---
    [pref, ~, Lref] = resample_closed_polyline(pref, opts.stepsize_opt);
    L = Lref;
    N = size(pref,1);

    % resample bounds to same N / station
    s_opt = (0:opts.stepsize_opt:(L - opts.stepsize_opt))';
    bL = interp_closed(bL, s_common, s_opt);
    bR = interp_closed(bR, s_common, s_opt);

    debug = struct(); debug.iters = [];
    alpha_prev = zeros(N,1);
    raceline_xy = pref;

    for it = 1:opts.max_iters

        % right-pointing normals (repo convention)
        [nx, ny] = normals_right_from_polyline(pref);

        % widths consistent with:
        %   boundR ≈ p + n*wR  => dot(boundR - p, n) > 0
        %   boundL ≈ p - n*wL  => dot(boundL - p, n) < 0
        [wL, wR] = widths_from_bounds_rightnormal(pref, [nx ny], bL, bR);

        % alpha bounds (paper eq. (3)): alpha in [-wL + wveh/2, wR - wveh/2]
        alpha_lb = -(wL - opts.width_opt/2);
        alpha_ub =  (wR - opts.width_opt/2);
        alpha_lb = min(alpha_lb, alpha_ub - 1e-6); % keep feasible

        % Build QP and solve
        qp = build_mincurv_qp(pref, [nx ny], alpha_lb, alpha_ub, opts);

        qp_opts = optimoptions('quadprog', ...
            'Algorithm','interior-point-convex', ...
            'Display', opts.quadprog_display);

        [alpha, ~, exitflag] = quadprog(qp.H, qp.f, qp.A, qp.b, [], [], qp.lb, qp.ub, alpha_prev, qp_opts);

        if exitflag <= 0
            warning('quadprog did not converge (exitflag=%d). Clipping previous alpha.', exitflag);
            alpha = max(min(alpha_prev, alpha_ub), alpha_lb);
        end

        % IQP early scaling (paper: 1/3 and 2/3 for iter 1 and 2)
        if it == 1
            alpha = opts.alpha_scale_iter1 * alpha;
        elseif it == 2
            alpha = opts.alpha_scale_iter2 * alpha;
        end

        % Update raceline
        raceline_xy = pref + [nx ny].*alpha;

        % keep equal step size (paper mentions spline interpolation to keep equal step sizes)
        if opts.do_resample_each_iter
            [raceline_xy, ~, ~] = resample_closed_polyline(raceline_xy, opts.stepsize_opt);
            if size(raceline_xy,1) ~= N
                s_tmp = (0:opts.stepsize_opt:(opts.stepsize_opt*(size(raceline_xy,1)-1)))';
                raceline_xy = interp_closed(raceline_xy, s_tmp, s_opt);
            end
        end

        % Termination criterion: curvature profiles computed from linearizations
        dk = curvature_linearization_gap(pref, raceline_xy, opts.stepsize_opt);
        debug.iters(end+1) = struct('iter',it,'max_dkappa',dk,'alpha_norm',norm(alpha)); %#ok<AGROW>

        if it >= opts.min_iters && dk <= opts.kappa_error_allowed
            break;
        end

        % replace reference line by solution (IQP)
        pref = raceline_xy;
        alpha_prev = alpha;
    end

    % output curvature
    kappa = curvature_from_periodic_spline(raceline_xy, opts.stepsize_opt);

    % output station and optional final resample
    if opts.stepsize_out ~= opts.stepsize_opt
        [raceline_xy, ~, Lout] = resample_closed_polyline(raceline_xy, opts.stepsize_out); %#ok<ASGLU>
        kappa = curvature_from_periodic_spline(raceline_xy, opts.stepsize_out);
        Nout = size(raceline_xy,1);
        s_out = (0:opts.stepsize_out:(opts.stepsize_out*(Nout-1)))';
    else
        s_out = s_opt;
    end

    % final reftrack in repo format [x y w_right w_left]
    [nx, ny] = normals_right_from_polyline(pref);
    [wL, wR] = widths_from_bounds_rightnormal(pref, [nx ny], bL, bR);
    reftrack = [pref(:,1) pref(:,2) wR wL];

    out = struct();
    out.reftrack    = reftrack;
    out.raceline_xy = raceline_xy;
    out.alpha       = alpha;
    out.kappa       = kappa;
    out.s           = s_out;
    out.debug       = debug;
end

% ========================= DEFAULT OPTIONS =========================
function opts = set_default_opts(opts)
    if ~isfield(opts,'stepsize_prep'), opts.stepsize_prep = 1.0; end
    if ~isfield(opts,'stepsize_opt'),  opts.stepsize_opt  = 3.0; end
    if ~isfield(opts,'stepsize_out'),  opts.stepsize_out  = 2.0; end
    if ~isfield(opts,'width_opt'),     opts.width_opt     = 3.4; end
    if ~isfield(opts,'curvlim'),       opts.curvlim       = 0.12; end
    if ~isfield(opts,'use_curv_constr'), opts.use_curv_constr = true; end

    if ~isfield(opts,'max_iters'), opts.max_iters = 10; end
    if ~isfield(opts,'min_iters'), opts.min_iters = 3; end

    % repo default is 0.01; paper example uses 0.005 on Berlin
    if ~isfield(opts,'kappa_error_allowed'), opts.kappa_error_allowed = 0.01; end

    if ~isfield(opts,'alpha_scale_iter1'), opts.alpha_scale_iter1 = 1/3; end
    if ~isfield(opts,'alpha_scale_iter2'), opts.alpha_scale_iter2 = 2/3; end

    if ~isfield(opts,'do_smooth'), opts.do_smooth = true; end
    if ~isfield(opts,'smooth_window'), opts.smooth_window = 21; end
    if ~isfield(opts,'do_resample_each_iter'), opts.do_resample_each_iter = true; end

    if ~isfield(opts,'reg_H'), opts.reg_H = 1e-9; end
    if ~isfield(opts,'quadprog_display'), opts.quadprog_display = 'off'; end
end

% ========================= GEOMETRY HELPERS =========================
function xy = ensure_closed(xy)
    if norm(xy(1,:) - xy(end,:)) > 1e-9
        xy = [xy; xy(1,:)];
    end
end

function [xy_rs, s, L] = resample_closed_polyline(xy, ds)
    xy = ensure_closed(xy);
    d = sqrt(sum(diff(xy,1,1).^2,2));
    s0 = [0; cumsum(d)];
    L  = s0(end);
    s = (0:ds:(L - ds))';
    xy_rs = interp_closed(xy, s0, s);
end

function xyq = interp_closed(xy_support, s_support, sq)
    xq = interp1(s_support, xy_support(:,1), sq, 'pchip');
    yq = interp1(s_support, xy_support(:,2), sq, 'pchip');
    xyq = [xq yq];
end

function xy_sm = smooth_closed(xy, win)
    if win < 3, xy_sm = xy; return; end
    if mod(win,2) == 0, win = win + 1; end
    xy_sm = [circular_movmean(xy(:,1),win), circular_movmean(xy(:,2),win)];
end

function ysm = circular_movmean(y, win)
    n = numel(y);
    k = floor(win/2);
    yext = [y(end-k+1:end); y; y(1:k)];
    ker = ones(win,1)/win;
    yflt = conv(yext, ker, 'same');
    ysm = yflt(k+1:k+n);
end

function [nx, ny] = normals_right_from_polyline(xy)
    % tangent by central differences (circular)
    n = size(xy,1);
    ip = mod((1:n), n) + 1;
    im = mod((0:n-1) - 1, n) + 1;

    dx = xy(ip,1) - xy(im,1);
    dy = xy(ip,2) - xy(im,2);

    tnorm = sqrt(dx.^2 + dy.^2) + 1e-12;
    tx = dx ./ tnorm;
    ty = dy ./ tnorm;

    % RIGHT normal (repo convention): [ty, -tx]
    nx =  ty;
    ny = -tx;
end

function [wL, wR] = widths_from_bounds_rightnormal(pref, nvec, bL, bR)
    nx = nvec(:,1); ny = nvec(:,2);
    vR = bR - pref;
    vL = bL - pref;

    projR = vR(:,1).*nx + vR(:,2).*ny; % should be >0
    projL = vL(:,1).*nx + vL(:,2).*ny; % should be <0

    % if signs look swapped, swap boundaries
    if mean(projR) < 0 && mean(projL) > 0
        tmp = bL; bL = bR; bR = tmp;
        vR = bR - pref; vL = bL - pref;
        projR = vR(:,1).*nx + vR(:,2).*ny;
        projL = vL(:,1).*nx + vL(:,2).*ny;
    end

    wR = max(projR, 0);
    wL = max(-projL, 0);
end

% ========================= SPLINE-BASED QP BUILD =========================
function qp = build_mincurv_qp(pref, nvec, alpha_lb, alpha_ub, opts)
    N = size(pref,1);
    px = pref(:,1); py = pref(:,2);
    nx = nvec(:,1); ny = nvec(:,2);

    h = ones(N,1) * opts.stepsize_opt;

    S = spline_second_derivative_operator_periodic(h);

    xdd0 = S * px;
    ydd0 = S * py;

    Sx = S * diag(nx);
    Sy = S * diag(ny);

    xprime = spline_first_derivative_at_knots(px, h, xdd0);
    yprime = spline_first_derivative_at_knots(py, h, ydd0);

    den3  = (xprime.^2 + yprime.^2).^3 + 1e-18;
    Pxx = spdiags((yprime.^2)./den3, 0, N, N);
    Pxy = spdiags((-2*xprime.*yprime)./den3, 0, N, N);
    Pyy = spdiags((xprime.^2)./den3, 0, N, N);

    Hq = (Sx' * Pxx * Sx) + (Sy' * Pxy * Sx) + (Sy' * Pyy * Sy);
    fq = (2*Sx' * Pxx * xdd0) + (Sy' * Pxy * xdd0) + (Sx' * Pxy * ydd0) + (2*Sy' * Pyy * ydd0);

    H = 2 * Hq;              % quadprog uses 1/2*x'*H*x
    H = 0.5*(H + H');
    H = H + opts.reg_H * speye(N);

    A = []; b = [];
    if opts.use_curv_constr
        den32 = (xprime.^2 + yprime.^2).^(3/2) + 1e-18;
        Qx = spdiags(yprime ./ den32, 0, N, N);
        Qy = spdiags(xprime ./ den32, 0, N, N);

        k_ref = (Qy*ydd0) - (Qx*xdd0);
        E     = (Qy*Sy)   - (Qx*Sx);

        A = [E; -E];
        b = [opts.curvlim - k_ref; opts.curvlim + k_ref];
    end

    qp = struct('H',H,'f',fq,'A',A,'b',b,'lb',alpha_lb,'ub',alpha_ub);
end

function S = spline_second_derivative_operator_periodic(h)
    N = numel(h);
    ip = mod((1:N), N) + 1;
    im = mod((0:N-1) - 1, N) + 1;

    A = spalloc(N, N, 3*N);
    B = spalloc(N, N, 3*N);

    for i = 1:N
        hm = h(im(i));
        hi = h(i);

        A(i, im(i)) = hm;
        A(i, i)     = 2*(hm + hi);
        A(i, ip(i)) = hi;

        B(i, im(i)) =  6/hm;
        B(i, i)     = -6*(1/hm + 1/hi);
        B(i, ip(i)) =  6/hi;
    end

    S = A \ B;
end

function d = spline_first_derivative_at_knots(y, h, m)
    N = numel(y);
    ip = mod((1:N), N) + 1;
    d = (y(ip) - y) ./ h - (h .* (2*m + m(ip)))/6;
end

% ========================= CURVATURE + TERMINATION =========================
function kappa = curvature_from_periodic_spline(xy, ds)
    N = size(xy,1);
    h = ones(N,1)*ds;
    S = spline_second_derivative_operator_periodic(h);

    x = xy(:,1); y = xy(:,2);
    xdd = S*x; ydd = S*y;

    xp = spline_first_derivative_at_knots(x, h, xdd);
    yp = spline_first_derivative_at_knots(y, h, ydd);

    den = (xp.^2 + yp.^2).^(3/2) + 1e-18;
    kappa = (xp.*ydd - yp.*xdd) ./ den;
end

function dk = curvature_linearization_gap(ref_xy, sol_xy, ds)
    % Compare curvature computed with "linearization along ref" vs "along sol"
    N = size(ref_xy,1);
    h = ones(N,1)*ds;
    S = spline_second_derivative_operator_periodic(h);

    % sol second derivatives
    xdd_s = S*sol_xy(:,1);
    ydd_s = S*sol_xy(:,2);

    % ref first derivatives -> Qref
    xdd_r = S*ref_xy(:,1);
    ydd_r = S*ref_xy(:,2);
    xp_r  = spline_first_derivative_at_knots(ref_xy(:,1), h, xdd_r);
    yp_r  = spline_first_derivative_at_knots(ref_xy(:,2), h, ydd_r);
    den_r = (xp_r.^2 + yp_r.^2).^(3/2) + 1e-18;
    Qx_r = spdiags(yp_r ./ den_r, 0, N, N);
    Qy_r = spdiags(xp_r ./ den_r, 0, N, N);

    k_refLin = (Qy_r*ydd_s) - (Qx_r*xdd_s);

    % sol first derivatives -> Qsol
    xp_s = spline_first_derivative_at_knots(sol_xy(:,1), h, xdd_s);
    yp_s = spline_first_derivative_at_knots(sol_xy(:,2), h, ydd_s);
    den_s = (xp_s.^2 + yp_s.^2).^(3/2) + 1e-18;
    Qx_s = spdiags(yp_s ./ den_s, 0, N, N);
    Qy_s = spdiags(xp_s ./ den_s, 0, N, N);

    k_solLin = (Qy_s*ydd_s) - (Qx_s*xdd_s);

    dk = max(abs(k_solLin - k_refLin));
end