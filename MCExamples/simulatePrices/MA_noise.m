function [chi] = MA_noise(num,rep,d,M,sigma_0)
% Reference: Eq. 4.15 in Jacod, Li & Zheng (2019, Journal of Econometrics)

chi = zeros(num,rep);
Z = sigma_0*randn(num,rep);
phi = zeros(M,1);

for i = 1:M
    temp = (d:1:i-1+d);
    phi(i) = prod(temp)/factorial(i);
end

chi(1,:) = Z(1,:);
for i = 2:num
    temp = zeros(1,rep);
    if i <= M
        for m = 1:i-1
            temp = temp + phi(m)*Z(i-m,:);
        end
    else
        for m = 1:M
            temp = temp + phi(m)*Z(i-m,:);
        end
    end
    chi(i,:) = Z(i,:) - temp;
end

end
