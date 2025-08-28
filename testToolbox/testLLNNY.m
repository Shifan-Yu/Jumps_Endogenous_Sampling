function [test, m_hat_eps, m_hat, N_c] = testLLNNY(price, c, epsilon, H2_tab, H2eps_tab, AVAR)
% testLLNNY — Test statistic with implied barriers from h_2(m) and h_2,eps(m)
% -----------------------------------------------------------------------------------
%   • Path X by 'price'.
%   • Barrier c determines PDS hits and returns r^{(c)}.
%   • h_2(m)   = E[(r^{(c)})^2] / c^2.
%   • h_2,ε(m) = E[ min{ (r^{(c)})^2, (c(1+ε))^2 } ] / c^2.
%   • Tables provided on an m-grid (two columns): 
%       H2_tab   = [m_grid, h_2(m_grid)]
%       Heps_tab = [m_grid, h^ε(m_grid)]     (for the SAME ε)
%       AVAR     = [m_grid, AVAR(m_grid)]    (Asymptotic variance for the SAME ε)
%
% Steps:
%   1) Extract PDS returns at c: r^{(c)} and compute squared returns r2.
%   2) Empirical scale-free moments:
%         Y_0      = (1/N_c) Σ r2 / c^2                     ≈ h_2(m)
%         Y_eps    = (1/N_c) Σ min(r2, (c(1+ε))^2) / c^2    ≈ h^ε(m)
%   3) Invert the monotone maps on the grid:
%         m_hat     : h_2(m_hat)     ≈ Y_0
%         m_hat_eps : h^ε(m_hat_eps) ≈ Y_eps
%   4) Standardize T_LLNNY = ((m_hat_eps/m_hat) − 1) by √( AVAR(m_hat_eps) / n )
% -----------------------------------------------------------------------------------

    % 1) PDS returns at barrier m, then square
    r2 = ret_delta(price, c).^2;

    % N_m = number of PDS hits at barrier m
    N_c = length(r2);
    if N_c == 0
        test = NaN; m_hat_eps = NaN; m_hat = NaN;
        warning('No PDS hits at c=%.6g. Consider a smaller barrier.', c);
        return;
    end

    % 2) Empirical scale-free moments (match definitions of h_2 and h^ε)
    cap_level = (c * (1 + epsilon))^2;      % (m(1+ε))^2
    r2_cens   = min(r2, cap_level);

    Y_0   = sum(r2)      / (c^2 * N_c);     % ≈ h_2(m)
    Y_eps = sum(r2_cens) / (c^2 * N_c);     % ≈ h^ε(m)

    % 3) Implied barriers via monotone inversion on the tabulated grid
    m_hat     = invFunc(Y_0,   H2_tab);     % h_2(m_hat)     ≈ Y_0
    m_hat_eps = invFunc(Y_eps, H2eps_tab);  % h^ε(m_hat_eps) ≈ Y_eps

    % 4) AVAR evaluated at m_hat_eps via linear interpolation
    AVAR_m = linearInt(m_hat_eps, AVAR);

    % LLNNY statistic
    n = length(price);
    test = (m_hat_eps / m_hat - 1) / sqrt( max(realmin, AVAR_m) / n );

end

