%% Simulate mu2(x), h(x)=mu2/x^2, epsilon-censored moments, derivatives, and AVAR
% -----------------------------------------------------------------------------------------
% This script reproduces (via Monte Carlo) the paper’s objects:
%   • PDS return r^{(x)} at symmetric barrier x for a GRW (Gaussian random walk).
%   • mu2(x)   = E[(r^{(x)})^2]                                  — 2nd moment of PDS returns.
%   • h(x)     = mu2(x) / x^2  (= h_2(m) in the paper)
%   • Censoring: 
%       r_c^{(x,eps)} = min{ (r^{(x)})^2, (x(1+eps))^2 }
%       mu2_eps(x,eps) = E[r_c^{(x,eps)}]
%       h_eps(x,eps) = mu2_eps / x^2  (= h_2,eps(m) in the paper)
%   • Derivatives: h'(x), h'_eps(x,eps) via local linear OLS slope.
%   • AVAR formulas (for linear and ratio tests, respectively) from the delta-method expressions in the paper:
%       AVAR_linear(x,eps)  = mu2(x) * [ v(x)/h'(x)^2 + v_eps(x,eps)/h'_eps(x,eps)^2
%                                      - 2*c_eps(x,eps)/(h'(x) h'_eps(x,eps)) ] / x^4
%       AVAR_ratio(x,eps) = AVAR_left(x,eps) / x^2
% -----------------------------------------------------------------------------------------

clearvars; clc;
% parpool(16);

rng(1);  % Reproducibility

%% 1) Simulate Gaussian Random Walk (GRW)
% For very large x you’ll need a much larger 'num' (and likely ret_delta_mex).
num = 1e6;                                  % demo run (fast)
% num = 1e9;                                % what we actually use! 
GRW = cumsum(randn(num,1));

%% 2) Grids and preallocations (keep original variable names)
x_all   = (0.01:0.01:90)';                  % barrier grid for PDS (paper: m-grid)
eps_all = (0.01:0.01:1)';                   % censoring levels (0,1]

mu2 = zeros(length(x_all),2); mu2(:,1) = x_all;      % [x, mu2(x)]
v   = zeros(length(x_all),2); v(:,1)   = x_all;      % [x, Var((r^{(x)})^2)]
h   = zeros(length(x_all),2); h(:,1)   = x_all;      % [x, h(x)=mu2/x^2]

mu2_eps = zeros(length(x_all),length(eps_all)+1); mu2_eps(:,1) = x_all; % [x, mu2_eps(x,eps)]
v_eps   = zeros(length(x_all),length(eps_all)+1);   v_eps(:,1) = x_all; % [x, Var(censored)]
c_eps   = zeros(length(x_all),length(eps_all)+1);   c_eps(:,1) = x_all; % [x, Cov((r^2), censored)]
h_eps   = zeros(length(x_all),length(eps_all)+1);   h_eps(:,1) = x_all; % [x, h_eps(x,eps)]

%% 3) Core loop over x: (r^{(x)})^2, epsilon-censored moments, and h(x)
% ret_delta(GRW,x) must return the PDS return series for barrier x.
% We square immediately to align with mu2 and AVAR formulas in the paper.
E = (1 + eps_all').^2;  % 1-by-Neps; cap level for (r^{(x)})^2 at (x^2)*(1+eps)^2.

for i = 1:length(x_all)
    x = x_all(i);

    % PDS returns at barrier x -> squared returns drive mu2(x)
    % Swap to ret_delta_mex(GRW,x).^2 if you have the MEX version.
    ret = ret_delta(GRW, x).^2;

    % Edge case: too-large x may produce zero samples for short paths
    if isempty(ret)
        warning('Empty PDS sample at x=%.4f. Increase ''num'' or reduce x.', x);
        continue;
    end

    % Uncensored moment/variance and scale-free map
    mu2(i,2) = mean(ret);
    v(i,2)   = var(ret);            % unbiased sample variance
    h(i,2)   = mu2(i,2) / x^2;      % h_2(m) in the paper

    % Vectorized epsilon-censoring across all eps (robust to MATLAB versions)
    % rc(j,k) = min( ret(j), (x^2)*E(k) )  =>  size(ret)=n x 1, E=1 x Neps
    rc = bsxfun(@min, ret, (x^2) * E);   % length(ret)-by-Neps matrix

    % First/second moments of censored variable
    mu2_eps(i,2:end) = mean(rc, 1);
    v_eps(i,2:end)   = var(rc, 0, 1);    % unbiased along dim=1

    % Sample covariance Cov( (r^{(x)})^2, rc )
    rc_centered      = bsxfun(@minus, rc, mu2_eps(i,2:end));
    c_eps(i,2:end)   = (ret - mean(ret))' * rc_centered / (length(ret) - 1);

    % Scale-free censored map
    h_eps(i,2:end)   = mu2_eps(i,2:end) / x^2;

    % disp(x); % optional progress
end

% % Optional: persist base tables
% save('mu2_vec.mat','mu2');
% save('mu2_eps_vec.mat','mu2_eps');
% save('v_vec.mat','v');
% save('v_eps_vec.mat','v_eps');
% save('c_eps_vec.mat','c_eps');
% save('h_vec.mat','h');
% save('h_eps_vec.mat','h_eps');

%% 4) Derivatives of h(x) and h_eps(x,eps): local linear slope
% We keep evaluation inside the x_all grid to avoid extrapolation and
% boundary artifacts (which inflate noise in the slope estimator).
x_lo   = max(x_all(1),  0.01);
x_hi   = min(x_all(end) - 0.01, 89);   % slightly interior to avoid edges
xrange = (0.01:0.01:x_hi)';              % keep your left bound; clamp the right

% Align rows by mapping xrange to indices in x_all (protect against fp drift)
[~, idx] = ismember(round(xrange, 5), round(x_all, 5));
if ~all(idx > 0)
    error('xrange must be a subset of x_all; adjust the bounds.');
end

h_dev     = zeros(length(xrange),2);                  h_dev(:,1)     = xrange;  % [x, h'(x)]
h_eps_dev = zeros(length(xrange),length(eps_all)+1);  h_eps_dev(:,1) = xrange;  % [x, h'_eps(x,eps)]

% h'(x) via local-linear OLS on a ~100-point window (helper function)
for i = 1:length(xrange)
    xtrial = xrange(i);
    h_dev(i,2) = h_first_derivative(xtrial, h);              % uses [x, h(x)]
end

% h'_eps(x,eps_k) via same method, applied per epsilon column
for k = 1:length(eps_all)
    hk = h_eps(:,[1, k+1]);                                  % [x, h_eps_k(x)]
    for i = 1:length(xrange)
        xtrial = xrange(i);
        h_eps_dev(i,k+1) = h_first_derivative(xtrial, hk);
    end
end

% save('h_dev.mat','h_dev');
% save('h_eps_dev.mat','h_eps_dev');

%% 5) AVAR (left/right) on aligned xrange (delta-method expressions)
m = length(xrange);
mu2_x   = mu2(idx,2);
v_x     = v(idx,2);
v_eps_x = v_eps(idx,:);
c_eps_x = c_eps(idx,:);

avar_l = zeros(m, length(eps_all)+1); avar_l(:,1) = xrange;
avar_r = zeros(m, length(eps_all)+1); avar_r(:,1) = xrange;

for k = 2:length(eps_all)+1
    denom_h    = h_dev(:,2);      % h'(x)
    denom_heps = h_eps_dev(:,k);  % h'_eps(x,eps_k)

    % Guard: near-zero derivatives can explode AVAR (flat h around large x)
    small = (abs(denom_h) < 1e-12) | (abs(denom_heps) < 1e-12);
    if any(small)
        warning('Near-zero derivative at some x; AVAR may be unstable.');
    end

    avar_l(:,k) = mu2_x .* ( ...
        v_x ./ (denom_h.^2) + ...
        v_eps_x(:,k) ./ (denom_heps.^2) - ...
        2*c_eps_x(:,k) ./ (denom_h .* denom_heps) ) ./ (xrange.^4);

end

avar_r(:,2:end) = avar_l(:,2:end) ./ (xrange.^2);

% save('avar_l.mat','avar_l');
% save('avar_r.mat','avar_r');