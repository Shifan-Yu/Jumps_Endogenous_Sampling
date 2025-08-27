function res = linearInt(x_trial, xy_vec, method)
% linearInt  Linear interpolation y(x_trial) from a two-column table [x, y].
% Wraps MATLAB interp1 with safe sorting, NaN cleaning, and edge clamping.
% ------------------------------------------------------------------------------
% INPUTS
%   x_trial : scalar query point
%   xy_vec  : [N x 2] table with x in col 1 and y in col 2
%   method  : (optional) 'linear' (default) or 'pchip'
%
% OUTPUT
%   res     : interpolated value at x_trial (clamped to ends if out-of-range)
% ------------------------------------------------------------------------------

    if nargin < 3 || isempty(method), method = 'linear'; end

    % Split and clean
    x = xy_vec(:,1);
    y = xy_vec(:,2);
    good = isfinite(x) & isfinite(y);
    x = x(good); y = y(good);

    if isempty(x)
        res = NaN; 
        return;
    end

    % Sort by x ascending (interp1 expects monotone x)
    [x, ix] = sort(x, 'ascend');
    y = y(ix);

    % Clamp the query to avoid extrapolation surprises
    xq = min(max(x_trial, x(1)), x(end));

    % Interpolate
    res = interp1(x, y, xq, method);   % 'linear' matches your original helper behavior
end
