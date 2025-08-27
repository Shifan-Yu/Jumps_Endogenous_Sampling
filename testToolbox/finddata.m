function [H2eps_tab, AVAR_tab] = finddata(eps, h_eps, AVAR)
% finddata  Extract the two-column tables for a specific ε from wide matrices.
% Assumes ε ∈ {0.01,0.02,...,1.00} and column k+1 holds ε=k/100.
    col = max(2, min(size(h_eps,2), floor(eps*100) + 1));
    H2eps_tab = h_eps(:, [1, col]);       % [m, h_2,ε(m)]
    AVAR_tab  = AVAR(:, [1, col]);      % [m, AVAR(m,ε)]
end