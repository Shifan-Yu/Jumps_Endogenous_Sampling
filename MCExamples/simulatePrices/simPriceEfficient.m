function [x,xc,xd,v,N,iv] = simPriceEfficient...
    (num,rep,initial,mu,kappa,theta,xi,rho,lambda,beta)

dt = zeros(num,rep) + 1/num;
xc = zeros(num+1,rep); 
xd = zeros(num+1,rep); 
v = zeros(num+1,rep);
xc(1,:) = initial;
v(1,:) = theta;
dW1 = randn(num,rep);
dW = randn(num,rep);
dW2 = dW1*rho + dW*sqrt(1-rho^2);

for i = 1:num
    v(i+1,:) = v(i,:) + kappa*(theta-v(i,:)).*dt(i,:) + xi*sqrt(max(0,v(i,:))).*sqrt(dt(i,:)).*dW1(i,:);
    xc(i+1,:) = xc(i,:) + mu.*dt(i,:) + sqrt(max(0,v(i,:))).*sqrt(dt(i,:)).*dW2(i,:);
end

N = zeros(rep,1);
parfor m = 1:rep
    xdtemp = zeros(num+1,1);
    Jloc = find(rand([num,1])<lambda/num); 
    if ~isempty(Jloc)
        N(m) = length(Jloc);
        for i = 1:length(Jloc)
            loc = Jloc(i)+1;
            Jsize = beta*randLaplace(1,1)
            % Laplace distribution (centred) with the scale parameter beta
            xdtemp(loc:num+1) = xdtemp(loc:num+1) + sign(rand(1)-0.5)*Jsize;
        end
    end
    xd(:,m) = xdtemp;
end

x = xc + xd;
iv = sum(v/num)';
    
end