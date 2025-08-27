function [ret] = ret_delta(price,delta)

m = 1;
T = length(price);
ret = zeros(T,1);
for i = 2:T
    if abs(price(i)-price(m)) >= delta
        ret(m) = price(i)-price(m);
        m = i;
    end
end
ret(ret == 0) = [];

end
