function par = ena_param_struct_from_veh(veh)
%ena_PARAM_STRUCT_FROM_VEH Map vehicle parameters to vehicle9_ena_eom_autogen par struct.
%
% The ena EOM uses a sprung rigid body plus four rigidly mounted
% unsprung point masses.  mWf and mWr are per-wheel unsprung masses.

    m      = getvehfield_any_local(veh, {'m'}, 1450.0);
    m_uf   = getvehfield_any_local(veh, {'m_uf','m_unsprung_f','muf','m_uns_f'}, 0.0);
    m_ur   = getvehfield_any_local(veh, {'m_ur','m_unsprung_r','mur','m_uns_r'}, 0.0);
    m_s    = getvehfield_any_local(veh, {'m_s','m_sprung','ms'}, m - m_uf - m_ur);

    acom   = getvehfield_any_local(veh, {'acom','lf'}, 1.35);
    bcom   = getvehfield_any_local(veh, {'bcom','lr'}, 1.35);
    wt     = getvehfield_any_local(veh, {'wt'}, 1.60);
    tF     = getvehfield_any_local(veh, {'wt_f','tF','track_f','tf'}, wt);
    tR     = getvehfield_any_local(veh, {'wt_r','tR','track_r','tr'}, wt);

    Iz_total = getvehfield_any_local(veh, {'Iz'}, 2200.0);
    unsprung_yaw = m_uf*(acom^2 + (tF/2)^2) + m_ur*(bcom^2 + (tR/2)^2);

    par = struct();
    par.m_s = m_s;
    par.mWf = 0.5*m_uf;
    par.mWr = 0.5*m_ur;
    par.Ixx = getvehfield_any_local(veh, {'Ix_sprung','Ix_s','Ixs','Ix'}, 650.0);
    par.Iyy = getvehfield_any_local(veh, {'Iy_sprung','Iy_s','Iys','Iy'}, 1500.0);
    par.Izz = getvehfield_any_local(veh, {'Izz_sprung','Iz_sprung','Izz_s'}, Iz_total - unsprung_yaw);
    par.Ixz = getvehfield_any_local(veh, {'Ixz'}, 0.0);
    par.Iw  = getvehfield_any_local(veh, {'Iw'}, 1.25);

    par.lf = acom;
    par.lr = bcom;
    par.tf = tF;
    par.tr = tR;
    par.hc = getvehfield_any_local(veh, {'hc','hcom'}, 0.55);
    par.rw = getvehfield_any_local(veh, {'rw'}, 0.31);

    par.Ktheta = getvehfield_any_local(veh, {'Ktheta'}, 85000.0);
    par.Ctheta = getvehfield_any_local(veh, {'Ctheta'}, 9000.0);
    par.Kphi   = getvehfield_any_local(veh, {'Kphi'}, 65000.0);
    par.Cphi   = getvehfield_any_local(veh, {'Cphi'}, 7500.0);

    if abs(m - (m_s + m_uf + m_ur)) > 1e-8
        warning('Mass split mismatch: m=%.12g but m_s+m_uf+m_ur=%.12g.', m, m_s + m_uf + m_ur);
    end
    if par.Izz <= 0
        warning('Computed ena sprung yaw inertia Izz=%.12g is non-positive. Check Iz and unsprung mass/geometry.', par.Izz);
    end
end

function val = getvehfield_any_local(veh, names, default_val)
    val = default_val;
    for i = 1:numel(names)
        if isfield(veh, names{i}) && ~isempty(veh.(names{i}))
            val = veh.(names{i});
            return;
        end
    end
end
