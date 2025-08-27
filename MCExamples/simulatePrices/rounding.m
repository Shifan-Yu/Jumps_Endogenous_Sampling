function [res] = rounding(price_vec)

res = log(round(exp(price_vec),2));

end
