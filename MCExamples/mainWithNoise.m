clearvars; clc;
% parpool(16);

%% ------------------------------------------------------------------------
%  DGP: Heston diffusion + compound Laplace jumps
%  X = X^c + X^d with:
%    dX^c_t =  mu dt + sqrt(v_t) dW^x_t,  
%    dv_t = kappa(theta - v_t)dt + xi sqrt(v_t) dW^v_t,
%    Corr(dW^x, dW^v) = rho,
%    X^d is stepwise with Bernoulli(λ/num) arrivals and Laplace(0,β) jump sizes.
% -------------------------------------------------------------------------
rep         = 500;             % # Monte Carlo paths (columns)
num         = 23400;           % points over [0,1]
initial     = 4;

% Heston params (dailyized)
mu          = 0.05/252;
kappa       = 5/252;
xi          = 0.5/252;
theta       = 0.16/252;
rho         = -0.5;

% Jumps: expected count λ, Laplace scale β
lambda      = 1;                        % E[# jumps over [0,1]]
beta1       = (1/5)*sqrt(theta);        % mild jumps
beta2       = 0.4*sqrt(theta);          % larger jumps

d           = 0.3;                      % MA_noise parameter
M           = 100;                      % MA_noise parameter
r           = 1;                        % rounding flag

% Simulate baseline diffusion + two jump scenarios (independent jump draws)
[x1, xc, xd1, v, N1, ~] = simPriceNoise_autoGau_t(num,rep,initial,mu,kappa,theta,xi,rho,lambda,beta1,d,M,r);
[~,  ~,  xd2, ~, N2, ~] = simPriceNoise_autoGau_t(num,rep,initial,mu,kappa,theta,xi,rho,lambda,beta2,d,M,r);
x2 = xc + xd2;                          % larger-jump alternative

%% ------------------------------------------------------------------------
%  Barrier family: c_i = K * sqrt(Var(ΔX_i))
% -------------------------------------------------------------------------

K_all = (3:1:10)';               % grid of barrier scalings K
eps   = 0.05;                    % ε used in the censoring threshold
kn    = round(0.5*sqrt(num));    % pre-averaging block size

load('h_vec.mat');               % h_vec      = [m_grid, h_2(m)]
load('h_eps_vec.mat');           % h_eps_vec  = [m_grid, h_2,ε(m)] (wide in ε)
load('avar_r.mat');              % avar_r     = [m_grid, AVAR_ε(m)] (wide in ε)
[H2eps_tab, AVAR_tab] = finddata(eps, h_eps_vec, avar_r);  % -> [m, h_2,ε(m)], [m, AVAR_ε(m)]
H2_tab = h_vec;                                            % -> [m, h_2(m)]

% Asymptotic (right-tail) critical value for size: Z ~ N(0,1)
zcrit_right = sqrt(2) * erfinv(2*0.95 - 1);

%% ------------------------------------------------------------------------
%  SIZE under the null (diffusion only: X = xc)
% -------------------------------------------------------------------------
Z_null  = NaN(rep, numel(K_all));      % test stats under null
Nc_null = zeros(rep, numel(K_all));    % # PDS hits (diagnostic)

for m = 1:numel(K_all)
    K = K_all(m);
    % parfor to speed up if a pool is open
    parfor i = 1:rep
        xi_path = wb_preaveraging(xc(:,i), kn); % pre-averaged path
        dX = diff(xi_path);
        c_i = K * sqrt(var(dX, 1));
        [Zi, ~, ~, Ni] = testLLNNY(xi_path, c_i, eps, H2_tab, H2eps_tab, AVAR_tab);
        Z_null(i,m)  = Zi;
        Nc_null(i,m) = Ni;
    end
end

% Rejection rate at 5% (omit NaNs)
size_rate = mean(Z_null > zcrit_right, 'omitnan');

% Empirical 95% and 99% null quantiles (per K) for size-adjusted power
q95 = quantile(Z_null, 0.95, 1);   

%% ------------------------------------------------------------------------
%  SIZE-ADJUSTED POWER (MILD JUMPS): keep only paths with N1>0
% -------------------------------------------------------------------------
idxA = find(N1 > 0);                   % keep only jumped paths
repA = numel(idxA);
Z_A  = NaN(repA, numel(K_all));

for m = 1:numel(K_all)
    K = K_all(m);
    parfor j = 1:repA
        i  = idxA(j);
        xi_path = wb_preaveraging(x1(:,i), kn); % pre-averaged path
        dX = diff(xi_path);
        c_i = K * sqrt(var(dX, 1));

        [Zi, ~, ~, ~] = testLLNNY(xi_path, c_i, eps, H2_tab, H2eps_tab, AVAR_tab);
        Z_A(j,m) = Zi;
    end
end

power_mild = mean(Z_A > q95, 'omitnan');   % compare to null’s 95% per K

%% ------------------------------------------------------------------------
%  SIZE-ADJUSTED POWER (LARGE JUMPS): keep only paths with N2>0
% -------------------------------------------------------------------------
idxB = find(N2 > 0);
repB = numel(idxB);
Z_B  = NaN(repB, numel(K_all));

for m = 1:numel(K_all)
    K = K_all(m);
    parfor j = 1:repB
        i  = idxB(j);
        xi_path = wb_preaveraging(x2(:,i), kn); % pre-averaged path
        dX = diff(xi_path);
        c_i = K * sqrt(var(dX, 1));

        [Zi, ~, ~, ~] = testLLNNY(xi_path, c_i, eps, H2_tab, H2eps_tab, AVAR_tab);
        Z_B(j,m) = Zi;
    end
end

power_large = mean(Z_B > q95, 'omitnan');

%% ------------------------------------------------------------------------
%  REPORT: per K, show size and size-adjusted power
% -------------------------------------------------------------------------
fprintf('\nε = %.3f, rep = %d\n', eps, rep);
fprintf('K    Size(5%%)   Power_mild   Power_large   N_c (null)\n');
for m = 1:numel(K_all)
    K = K_all(m);
    Nc = mean(Nc_null(:,m), 'omitnan');
    fprintf('%4.0f   %8.4f   %11.4f   %12.4f   %10.0f\n', ...
        K, size_rate(m), power_mild(m), power_large(m), Nc);
end
