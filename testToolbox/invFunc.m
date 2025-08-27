function m_hat = invFunc(Y_trial, mY_table)
% invFunc  Robust monotone inverse for tabulated Y(m) (e.g., h_2 or h^ε).
% ------------------------------------------------------------------------------
% Given Y_trial and a two-column table mY_table=[m_grid, Y_grid] with Y ↓ in m,
% return m_hat s.t. Y(m_hat) ≈ Y_trial.
% Two-stage approach:
%   (1) For deep-tail Y (very small Y / large m), fit μ2(m) on the tail
%       (since Y = μ2/m^2) via quadratic, then solve (a−Y_trial)m^2 + b m + c = 0.
%   (2) Else do a local inverse quadratic in u=1/m near the bracket of Y_trial.
% If roots are unusable, fall back to monotone PCHIP interpolation (Y→m).
% ------------------------------------------------------------------------------

    % --- Clean & extract inputs -------------------------------------------------
    m = mY_table(:,1);                    % m-grid (barrier)
    Y = mY_table(:,2);                    % Y(m) values
    good = isfinite(m) & isfinite(Y);     % drop NaNs
    m = m(good);  Y = Y(good);
    if numel(m) < 3, m_hat = NaN; return; end

    % --- Sort by m asc, dedupe, enforce Y decreasing ---------------------------
    [m, ix] = sort(m, 'ascend'); Y = Y(ix);
    [m, iu] = unique(m, 'stable'); Y = Y(iu);
    if Y(1) < Y(end)               % if not decreasing, flip so Y(1) is largest
        m = flipud(m); Y = flipud(Y);
    end

    % --- Clamp out-of-range Y_trial to grid ends -------------------------------
    if Y_trial >= Y(1)             % ask for too-large Y → smallest m
        m_hat = m(1);  return;
    end
    if Y_trial <= Y(end)           % ask for too-small Y → largest m
        m_hat = m(end); return;
    end

    % --- Try a tail quadratic only for deep-tail Y -----------------------------
    M = numel(m);
    tail_start = max(M - 10, 2);   % last ~10 points
    if Y_trial <= Y(tail_start)
        mt   = m(tail_start:end);
        mu2t = Y(tail_start:end) .* (mt.^2);      % since Y=μ2/m^2
        Xq   = [mt.^2, mt, ones(numel(mt),1)];
        abc  = Xq \ mu2t;                         % [a, b, c]
        A = (abc(1) - Y_trial);  B = abc(2);  C = abc(3);
        disc = B^2 - 4*A*C;
        if A ~= 0 && disc >= 0
            r1 = (-B + sqrt(disc)) / (2*A);
            r2 = (-B - sqrt(disc)) / (2*A);
            cand = [r1; r2];
            cand = cand(cand > 0 & isfinite(cand));
            if ~isempty(cand)
                % pick root closest to center of tail window, clamp to [min,max] of tail
                [~, j] = min(abs(cand - mean(mt)));
                m_hat  = min(max(cand(j), mt(1)), mt(end));
                return;                            % success via tail quad
            end
        end
        % else: fall through to local inverse
    end

    % --- Local inverse quadratic in u = 1/m near bracket of Y_trial -----------
    % Create increasing-Y copy to locate the bracket:
    Y_inc = flipud(Y);
    k = find(Y_inc >= Y_trial, 1, 'first');   % first index with Y_inc ≥ Y_trial
    center = M - k + 1;                       % map back to original indexing
    lo = max(1, center - 10); hi = min(M, center + 10);

    ml = m(lo:hi);     % local m window
    Yl = Y(lo:hi);     % local Y window
    u  = 1 ./ ml;      % regress in u=1/m (your original idea)
    Xl = [ones(numel(u),1), u, u.^2];

    b  = Xl \ Yl;                       % stable OLS: Y ≈ b0 + b1 u + b2 u^2
    b0 = b(1); b1 = b(2); b2 = b(3);

    % Solve b2 u^2 + b1 u + (b0 - Y_trial) = 0
    A2 = b2; B2 = b1; C2 = b0 - Y_trial;
    disc2 = B2^2 - 4*A2*C2;

    if abs(A2) < 1e-12                 % degenerate to linear in u
        if abs(B2) < 1e-12
            m_hat = interp1(flipud(Y), flipud(m), Y_trial, 'pchip'); % monotone fallback
            return;
        end
        u_hat = -C2 / B2;
    elseif disc2 >= 0
        u1 = (-B2 + sqrt(disc2)) / (2*A2);
        u2 = (-B2 - sqrt(disc2)) / (2*A2);
        uc = [u1; u2];                       % candidate u roots
        uc = uc(isfinite(uc) & uc > 0);      % keep positive finite
        if isempty(uc)
            m_hat = interp1(flipud(Y), flipud(m), Y_trial, 'pchip'); % fallback
            return;
        end
        mc = 1 ./ uc;                         % candidate m roots
        [~, j] = min(abs(mc - mean(ml)));     % pick near window center
        u_hat = uc(j);
    else
        m_hat = interp1(flipud(Y), flipud(m), Y_trial, 'pchip');     % complex roots → fallback
        return;
    end

    m_hat = 1 / u_hat;                        % back-transform
end
