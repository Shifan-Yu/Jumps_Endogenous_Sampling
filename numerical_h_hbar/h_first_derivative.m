function [h_dev] = h_first_derivative(x_trial, h_vec)
% h_first_derivative: local-linear slope estimate of h'(x) at x_trial
% ---------------------------------------------------------------------------------
% Method:
%   Fit OLS of h on x with an intercept over a local window around x_trial,
%   then take the slope as the derivative estimate.
%   This stabilizes the numerical derivative of h(x)=mu2(x)/x^2 for MC data.
%
% Inputs
%   x_trial : scalar, evaluation point
%   h_vec   : two-column array [x_all, h(x_all)]
%
% Output
%   h_dev   : estimated derivative h'(x_trial)
%
% ---------------------------------------------------------------------------------

    x_all = h_vec(:,1);
    h     = h_vec(:,2);
    n     = length(x_all);

    % Index of the rightmost grid point <= x_trial
    x_ind = sum(x_all <= x_trial);

    % Build a ~100-point window around x_ind, with caps at [1, n]
    if x_ind <= 0
        lo = 1;            hi = min(n, 100);
    elseif x_ind >= n
        lo = max(1, n-99); hi = n;
    else
        lo = max(1, x_ind-49);
        hi = min(n, x_ind+50);
    end
    rrange = lo:hi;

    % If the window degenerates (e.g., ultra-short grid), fall back gracefully
    if numel(rrange) < 2
        if x_ind >= 1 && x_ind < n
            dx   = x_all(x_ind+1) - x_all(x_ind);
            if dx == 0
                h_dev = 0;
            else
                h_dev = (h(x_ind+1) - h(x_ind)) / dx;
            end
        else
            h_dev = NaN;
        end
        return;
    end

    % Local linear regression (with intercept) over the window
    x_temp = [ones(numel(rrange),1), x_all(rrange)];
    h_temp = h(rrange);

    % Numerically stable OLS (avoids inv(X'X))
    beta_hat = x_temp \ h_temp;

    % Slope is the derivative estimate
    h_dev = beta_hat(2);
end

