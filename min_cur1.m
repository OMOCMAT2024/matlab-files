function out = mincurv_raceline_from_bounds(boundL_xy, boundR_xy, opts)
%MINCURV_RACELINE_FROM_BOUNDS
% Minimum-curvature raceline (QP + optional iterative re-linearization),
% based on Heilmeier et al. "Minimum curvature trajectory planning..." (VSD 2019)
% and the TUMFTM global_racetrajectory_optimization parameter conventions.
%
% INPUTS
%   boundL_xy : [NLx2] left boundary points [x y], ordered around the track (closed or open)
%   boundR_xy : [NRx2] right boundary points [x y], ordered around the track (closed or open)
%   opts      : struct with optional fields (defaults below)
%
% OUTPUT (struct "out")
%   out.reftrack      : [N x 4] [x y w_right w_left] reference line + widths (per station)
%   out.raceline_xy   : [N x 2] optimized raceline coordinates
%   out.alpha         : [N x 1] lateral shifts along normals (reference -> raceline)
%   out.kappa         : [N x 1] curvature along raceline (computed from periodic spline)
%   out.s             : [N x 1] station (0...L)
%   out.debug         : struct with extra iteration info
%
% REQUIREMENTS
%   Optimization Toolbox (quadprog)

    if nargin < 3, opts = struct(); end
    opts = set_default_opts(opts);

    % --- 1) Preprocess boundaries: close + resample to common station grid ---
    boundL_xy = ensure_closed(boundL_xy);
    boundR_xy = ensure_closed(boundR_xy);

    [bL_prep, sL, LL] = resample_closed_polyline(boundL_xy, opts.stepsize_prep);
    [bR_prep, sR, LR] = resample_closed_polyline(boundR_xy, opts.stepsize_prep);

    % use common length (assumes both boundaries represent the same loop)
    L = min(LL, LR);
    s_common = (0:opts.stepsize_prep:(L - opts.stepsize_prep))';
    bL = interp_closed(bL_prep, sL, s_common);
    bR = interp_closed(bR_prep, sR, s_common);

    % --- 2) Build initial reference line (midline) + smooth (optional) ---
    pref = 0.5 * (bL + bR); % initial reference points p_i
    if opts.do_smooth
        pref = smooth_closed(pref, opts.smooth_window);
    end

    % --- 3) Resample reference line to optimization stepsize (paper uses ~3m) ---
    [pref, sref, Lref] = resample_closed_polyline(pref, opts.stepsize_opt);
    L = Lref;
    N = size(pref,1);
    s = (0:opts.stepsize_opt:(L - opts.stepsize_opt))'; %#ok<NASGU> % length N

    % Also resample boundaries to same station count (for width initialization)
    bL = interp_closed(bL, s_common, (0:opts.stepsize_opt:(L - opts.stepsize_opt))');
    bR = interp_closed(bR, s_common, (0:opts.stepsize_opt:(L - opts.stepsize_opt))');

    % --- 4) Main iterative IQP loop (optional) ---
    debug = struct();
    debug.iters = [];

    alpha_prev = zeros(N,1);
    raceline_xy = pref;

    for it = 1:opts.max_iters

        % 4.1) periodic spline operators on current reference
        [nx, ny, tx, ty] = normals_from_polyline(pref); %#ok<ASGLU>

        % 4.2) widths to boundaries (projected along the (left) normal)
        [wL, wR, swapped] = widths_from_bounds(pref, [nx ny], bL, bR);
        if swapped
            tmp = bL; bL = bR; bR = tmp;
        end

        % bounds on alpha: alpha in [-w_left + wveh/2, w_right - wveh/2]
        % (paper Eq. (3))
        alpha_lb = -(wL - opts.width_opt/2);
        alpha_ub =  (wR - opts.width_opt/2);

        % clip widths if numeric issues
        alpha_lb = min(alpha_lb, alpha_ub - 1e-6);

        % 4.3) Build QP matrices from reference-line linearization
        qp = build_mincurv_qp(pref, [nx ny], alpha_lb, alpha_ub, opts);

        % 4.4) Solve QP
        qp_opts = optimoptions('quadprog', ...
            'Algorithm','interior-point-convex', ...
            'Display', opts.quadprog_display);

        [alpha, fval, exitflag] = quadprog(qp.H, qp.f, qp.A, qp.b, [], [], qp.lb, qp.ub, alpha_prev, qp_opts); %#ok<ASGLU>

        if exitflag <= 0
            warning('quadprog did not converge (exitflag=%d). Returning best available solution.', exitflag);
            alpha = max(min(alpha_prev, alpha_ub), alpha_lb);
        end

        % paper: scale alpha in iter 1 and 2 (1/3 and 2/3) to limit displacement
        if it == 1
            alpha = opts.alpha_scale_iter1 * alpha;
        elseif it == 2
            alpha = opts.alpha_scale_iter2 * alpha;
        end

        % 4.5) Update raceline and (optionally) re-sample to keep equal step size
        raceline_xy = pref + [nx ny].*alpha;
        if opts.do_resample_each_iter
            [raceline_xy, ~, ~] = resample_closed_polyline(raceline_xy, opts.stepsize_opt);
            % keep same N
            if size(raceline_xy,1) ~= N
                raceline_xy = interp_closed(raceline_xy, (0:opts.stepsize_opt:(opts.stepsize_opt*(size(raceline_xy,1)-1)))', ...
                                            (0:opts.stepsize_opt:(L - opts.stepsize_opt))');
            end
        end

        % 4.6) Termination criterion (paper: compare curvature linearizations)
        % "max difference between curvature profiles calculated based on linearizations
        %  along original reference line and along the result itself"
        [k_lin_ref, k_lin_self, k_true] = curvature_linearization_gap(pref, raceline_xy);

        dk = max(abs(k_lin_self - k_lin_ref));
        debug.iters(end+1) = struct('iter',it,'max_dkappa',dk,'alpha_norm',norm(alpha)); %#ok<AGROW>

        if it >= opts.min_iters && dk <= opts.kappa_error_allowed
            break;
        end

        % 4.7) Replace reference line by solution (paper iterative invocation)
        pref = raceline_xy;
        alpha_prev = alpha;
    end

    % final curvature of output raceline
    kappa = curvature_from_periodic_spline(raceline_xy);

    % resample final output if desired
    if opts.stepsize_out ~= opts.stepsize_opt
        [raceline_xy, ~, Lout] = resample_closed_polyline(raceline_xy, opts.stepsize_out); %#ok<ASGLU>
        kappa = curvature_from_periodic_spline(raceline_xy);
        Nout = size(raceline_xy,1);
        s_out = (0:opts.stepsize_out:(opts.stepsize_out*(Nout-1)))';
    else
        s_out = (0:opts.stepsize_opt:(L - opts.stepsize_opt))';
    end

    % build final reftrack format: [x y w_right w_left]
    % (repo uses [x y w_tr_right w_tr_left])
    [nx, ny] = normals_from_polyline(pref);
    [wL, wR] = widths_from_bounds(pref, [nx ny], bL, bR);
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
    % defaults inspired by TUMFTM racecar.ini (stepsize_prep=1m, stepsize_reg=3m, curvlim=0.12)
    if ~isfield(opts,'stepsize_prep'), opts.stepsize_prep = 1.0; end
    if ~isfield(opts,'stepsize_opt'),  opts.stepsize_opt  = 3.0; end
    if ~isfield(opts,'stepsize_out'),  opts.stepsize_out  = 2.0; end
    if ~isfield(opts,'width_opt'),     opts.width_opt     = 3.4; end  % vehicle width incl safety
    if ~isfield(opts,'curvlim'),       opts.curvlim       = 0.12; end % [rad/m]
    if ~isfield(opts,'use_curv_constr'), opts.use_curv_constr = true; end

    if ~isfield(opts,'max_iters'), opts.max_iters = 10; end
    if ~isfield(opts,'min_iters'), opts.min_iters = 3; end
    if ~isfield(opts,'kappa_error_allowed'), opts.kappa_error_allowed = 0.01; end

    if ~isfield(opts,'alpha_scale_iter1'), opts.alpha_scale_iter1 = 1/3; end
    if ~isfield(opts,'alpha_scale_iter2'), opts.alpha_scale_iter2 = 2/3; end

    if ~isfield(opts,'do_smooth'), opts.do_smooth = true; end
    if ~isfield(opts,'smooth_window'), opts.smooth_window = 21; end % odd integer recommended
    if ~isfield(opts,'do_resample_each_iter'), opts.do_resample_each_iter = true; end

    if ~isfield(opts,'reg_H'), opts.reg_H = 1e-9; end % small PSD regularization
    if ~isfield(opts,'quadprog_display'), opts.quadprog_display = 'off'; end
end

% ========================= GEOMETRY HELPERS =========================
function xy = ensure_closed(xy)
    if norm(xy(1,:) - xy(end,:)) > 1e-9
        xy = [xy; xy(1,:)];
    end
end

function [xy_rs, s, L] = resample_closed_polyline(xy, ds)
    % xy is closed (last == first). returns N points without duplicating last.
    xy = ensure_closed(xy);

    d = sqrt(sum(diff(xy,1,1).^2,2));
    s0 = [0; cumsum(d)];
    L  = s0(end);

    % remove the duplicated last point for interpolation, but keep it in "support"
    % support points for periodic interpolation:
    s_support = s0;
    xy_support = xy;

    % new station grid (exclude L)
    s = (0:ds:(L - ds))';

    xy_rs = interp_closed(xy_support, s_support, s);
end

function xyq = interp_closed(xy_support, s_support, sq)
    % xy_support includes last==first at s=L
    xq = interp1(s_support, xy_support(:,1), sq, 'pchip');
    yq = interp1(s_support, xy_support(:,2), sq, 'pchip');
    xyq = [xq yq];
end

function xy_sm = smooth_closed(xy, win)
    if win < 3
        xy_sm = xy; return;
    end
    if mod(win,2) == 0
        win = win + 1;
    end
    xy_sm = [circular_movmean(xy(:,1),win), circular_movmean(xy(:,2),win)];
end

function ysm = circular_movmean(y, win)
    % circular moving average for closed curve signals
    n = numel(y);
    k = floor(win/2);
    yext = [y(end-k+1:end); y; y(1:k)];
    ker = ones(win,1)/win;
    yflt = conv(yext, ker, 'same');
    ysm = yflt(k+1:k+n);
end

function [nx, ny, tx, ty] = normals_from_polyline(xy)
    % tangent by central differences (circular)
    n = size(xy,1);
    ip = mod((1:n), n) + 1;
    im = mod((0:n-1) - 1, n) + 1;

    dx = xy(ip,1) - xy(im,1);
    dy = xy(ip,2) - xy(im,2);

    tnorm = sqrt(dx.^2 + dy.^2) + 1e-12;
    tx = dx ./ tnorm;
    ty = dy ./ tnorm;

    % left normal
    nx = -ty;
    ny =  tx;
end

function [wL, wR, swapped] = widths_from_bounds(pref, nvec, bL, bR)
    % project boundary vectors onto current normal direction.
    % tries to auto-detect if L/R are swapped relative to driving direction.
    nx = nvec(:,1); ny = nvec(:,2);

    vL = bL - pref;
    vR = bR - pref;

    projL = vL(:,1).*nx + vL(:,2).*ny;   % should be positive (left)
    projR = vR(:,1).*nx + vR(:,2).*ny;   % should be negative (right)

    swapped = false;
    if mean(projL) < 0 && mean(projR) > 0
        % swap boundaries
        swapped = true;
        tmp = bL; bL = bR; bR = tmp;
        vL = bL - pref; vR = bR - pref;
        projL = vL(:,1).*nx + vL(:,2).*ny;
        projR = vR(:,1).*nx + vR(:,2).*ny;
    end

    wL = max(projL, 0);        % distance to left along +n
    wR = max(-projR, 0);       % distance to right along -n
end

% ========================= SPLINE-BASED QP BUILD =========================
function qp = build_mincurv_qp(pref, nvec, alpha_lb, alpha_ub, opts)
    % Build QP: minimize sum kappa^2 evaluated at knots, using periodic cubic spline.
    % Uses linear mapping: xdd = S * x, with x = px + nx .* alpha.
    % Then objective is quadratic in alpha.

    N = size(pref,1);
    px = pref(:,1); py = pref(:,2);
    nx = nvec(:,1); ny = nvec(:,2);

    % uniform spacing assumed (after resampling); still build h explicitly
    h = ones(N,1) * opts.stepsize_opt;

    % second-derivative operator for periodic cubic spline: m = S * y
    S = spline_second_derivative_operator_periodic(h);

    % reference second derivatives
    xdd0 = S * px;
    ydd0 = S * py;

    % alpha -> second derivatives
    Sx = S * diag(nx);
    Sy = S * diag(ny);

    % reference first derivatives (needed for P/Q)
    xprime = spline_first_derivative_at_knots(px, h, xdd0);
    yprime = spline_first_derivative_at_knots(py, h, ydd0);

    den3  = (xprime.^2 + yprime.^2).^3 + 1e-18;
    Pxx_d = (yprime.^2) ./ den3;
    Pxy_d = (-2*xprime.*yprime) ./ den3;
    Pyy_d = (xprime.^2) ./ den3;

    Pxx = spdiags(Pxx_d, 0, N, N);
    Pxy = spdiags(Pxy_d, 0, N, N);
    Pyy = spdiags(Pyy_d, 0, N, N);

    % quadratic term (cost = alpha' * Hq * alpha + fq' * alpha + const)
    Hq = (Sx' * Pxx * Sx) + (Sy' * Pxy * Sx) + (Sy' * Pyy * Sy);
    fq = (2*Sx' * Pxx * xdd0) + (Sy' * Pxy * xdd0) + (Sx' * Pxy * ydd0) + (2*Sy' * Pyy * ydd0);

    % quadprog wants 1/2*x'*H*x + f'*x
    H = 2 * (Hq);
    H = 0.5*(H + H'); % symmetrize
    H = H + opts.reg_H * eye(N);

    % constraints
    A = [];
    b = [];

    if opts.use_curv_constr
        % curvature kappa = (x' y'' - y' x'') / (x'^2+y'^2)^(3/2)
        den32 = (xprime.^2 + yprime.^2).^(3/2) + 1e-18;
        Qx = spdiags(yprime ./ den32, 0, N, N); % multiplies xdd with minus sign
        Qy = spdiags(xprime ./ den32, 0, N, N);

        k_ref = (Qy*ydd0) - (Qx*xdd0);
        E     = (Qy*Sy)   - (Qx*Sx);

        % |k_ref + E*alpha| <= curvlim  ->  E*alpha <= curvlim - k_ref  and  -E*alpha <= curvlim + k_ref
        A = [E; -E];
        b = [opts.curvlim - k_ref; opts.curvlim + k_ref];
    end

    qp = struct();
    qp.H  = H;
    qp.f  = fq;
    qp.A  = A;
    qp.b  = b;
    qp.lb = alpha_lb;
    qp.ub = alpha_ub;
end

function S = spline_second_derivative_operator_periodic(h)
    % Build S such that m = S*y gives periodic cubic spline second derivatives at knots.
    % Uses standard periodic spline system: A*m = B*y (cyclic tridiagonal).
    N = numel(h);
    ip = mod((1:N), N) + 1;
    im = mod((0:N-1) - 1, N) + 1;

    A = spalloc(N, N, 3*N);
    B = spalloc(N, N, 3*N);

    for i = 1:N
        hm = h(im(i));  % segment (i-1)->i
        hi = h(i);      % segment i->(i+1)

        A(i, im(i)) = hm;
        A(i, i)     = 2*(hm + hi);
        A(i, ip(i)) = hi;

        B(i, im(i)) =  6/hm;
        B(i, i)     = -6*(1/hm + 1/hi);
        B(i, ip(i)) =  6/hi;
    end

    S = full(A \ B);
end

function d = spline_first_derivative_at_knots(y, h, m)
    % derivative at start of each segment (knot i), periodic
    N = numel(y);
    ip = mod((1:N), N) + 1;

    d = (y(ip) - y) ./ h - (h .* (2*m + m(ip)))/6;
end

% ========================= CURVATURE + TERMINATION =========================
function kappa = curvature_from_periodic_spline(xy)
    % Compute curvature at knots from periodic cubic spline derivatives wrt station.
    N = size(xy,1);
    h = ones(N,1) * mean(sqrt(sum(diff([xy; xy(1,:)],1,1).^2,2))); %#ok<NASGU>
    % use uniform station spacing by re-resampling first for stability:
    % (here we assume xy is already uniform; if not, curvature quality degrades slightly)
    ds = mean(sqrt(sum(diff([xy; xy(1,:)],1,1).^2,2)));
    h = ones(N,1)*ds;

    S = spline_second_derivative_operator_periodic(h);

    x = xy(:,1); y = xy(:,2);
    xdd = S*x; ydd = S*y;

    xp = spline_first_derivative_at_knots(x, h, xdd);
    yp = spline_first_derivative_at_knots(y, h, ydd);

    den = (xp.^2 + yp.^2).^(3/2) + 1e-18;
    kappa = (xp.*ydd - yp.*xdd) ./ den;
end

function [k_lin_ref, k_lin_self, k_true] = curvature_linearization_gap(ref_xy, sol_xy)
    % k_lin_ref : curvature computed using (x',y') from ref, (x'',y'') from sol
    % k_lin_self: curvature computed using (x',y') from sol, (x'',y'') from sol
    % k_true    : same as k_lin_self here (true spline curvature at knots)

    N = size(ref_xy,1);
    ds = mean(sqrt(sum(diff([ref_xy; ref_xy(1,:)],1,1).^2,2)));
    h  = ones(N,1)*ds;

    S = spline_second_derivative_operator_periodic(h);

    % sol second derivatives
    xdd_s = S*sol_xy(:,1);
    ydd_s = S*sol_xy(:,2);

    % ref first derivatives
    xdd_r = S*ref_xy(:,1);
    ydd_r = S*ref_xy(:,2);
    xp_r  = spline_first_derivative_at_knots(ref_xy(:,1), h, xdd_r);
    yp_r  = spline_first_derivative_at_knots(ref_xy(:,2), h, ydd_r);

    den_r = (xp_r.^2 + yp_r.^2).^(3/2) + 1e-18;
    k_lin_ref = (xp_r.*ydd_s - yp_r.*xdd_s) ./ den_r;

    % sol first derivatives
    xp_s = spline_first_derivative_at_knots(sol_xy(:,1), h, xdd_s);
    yp_s = spline_first_derivative_at_knots(sol_xy(:,2), h, ydd_s);
    den_s = (xp_s.^2 + yp_s.^2).^(3/2) + 1e-18;

    k_lin_self = (xp_s.*ydd_s - yp_s.*xdd_s) ./ den_s;
    k_true = k_lin_self;
end