
function x = randLaplace(m,n)

u1 = rand(m,n);
u2 = rand(m,n);
x = log(u1./u2);

end