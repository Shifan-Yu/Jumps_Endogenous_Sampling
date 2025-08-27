function [x1] = wb_preaveraging(x,K)

[ret,~] = YbarYhat(diff(x),K,1);

rademacher = Rademacher_sim(length(ret));
ret = ret.*rademacher;
ret = ret(randperm(length(ret)));
x1 = cumsum(ret);

end

function [Ybar, Yhat] = YbarYhat(dY,kn,k)

% From https://sites.google.com/view/jiali/work

% [Ybar, Yhat] = YbarYhat(dY,kn,k)
% Compute Ybar, Yhat
%
% dY: the first difference of data
% kn: the window
% k:  the scale parameter

n = length(dY);
m = n - kn + 1;
Ybar = zeros(m,1);
Yhat = zeros(m,1);

% The weighting vector
g = gfun(k*(0:kn)/kn);
gprime2 = (g(2:end) - g(1:end-1)) .^ 2;
g(1) = [];

% Compute Ybar, Yhat
for i=1:m
    tmpbar = 0;
    tmphat = 0;
    for j=1:kn
        tmpbar = tmpbar + g(j) * dY(i+j-1);
        tmphat = tmphat + gprime2(j) * dY(i+j-1) ^ 2;
    end
    Ybar(i) = tmpbar;
    Yhat(i) = tmphat;
end
end

function output = gfun(x)
output = x .* (x>=0 & x<=0.5) + (1-x) .* (x<=1 & x>0.5);
end

function [x] = Rademacher_sim(n)

%n = number of samples; 
mp = [-1 1]; % <-- The two values you want as outputs
x = mp((rand(1,n)<.5)+1); % <-- Randomly pick one of the two values for n samples.
x = x';

end