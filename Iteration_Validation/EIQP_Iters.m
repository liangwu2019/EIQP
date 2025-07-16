function [Iters,Dualitys,Infeasibilitys] = EIQP_Iters(Q,c,A,b,epsilon)
nc = length(c);
nb = length(b);
n = nc+nb;

Iters = [];
Dualitys = [];
Infeasibilitys = [];

eta = 0.414213/sqrt(n+1);
gamma = 1-eta;

z = ones(nc,1);
y = ones(nb,1);
tau = 1;

v = ones(nc,1);
w = ones(nb,1);
kappa = 1;

Qz = Q*ones(nc,1);
ATy = A'*ones(nb,1);
Az = A*ones(nc,1);
zTc = ones(nc,1)'*c;
yTb = ones(nb,1)'*b;
zTQz = z'*Qz;

sigma = max([1;Qz-ATy+c; Az-b;-zTQz-zTc+yTb]);
Q = 1/sigma*Q;
A = 1/sigma*A;
c = 1/sigma*c;
b = 1/sigma*b; 

Qz = 1/sigma*Qz;
ATy = 1/sigma*ATy;
Az = 1/sigma*Az;
zTc = 1/sigma*zTc;
yTb = 1/sigma*yTb;
zTQz = 1/sigma*zTQz;

max_iter = ceil(log((n+1)/epsilon)/(-log(gamma)));

r_z = v - Qz + ATy -tau*c;
r_y = w - Az + tau*b;
r_tau = kappa + zTQz/tau + zTc - yTb;

M0 = [Q, -A', c;
      A, zeros(nb,nb), -b;
      -c', b', 0];

iter = 0;
while(1)
    iter = iter + 1;
    Iters = [Iters; iter];

    bar_mu = (z'*v+y'*w+tau*kappa)/(n+1);

    u_z = gamma*bar_mu./z - v + eta*r_z;
    u_y = gamma*bar_mu./y - w + eta*r_y;
    u_tau = gamma*bar_mu/tau - kappa + eta*r_tau;

    M = M0 + diag([v./z; w./y; zTQz/tau^2+kappa/tau]);
    M(n+1,1:nc) = M(n+1,1:nc) - 2*Qz'/tau;
    d_bar_x = M\([u_z;u_y;u_tau]);
    
    d_z = d_bar_x(1:nc);
    d_y = d_bar_x(nc+1:nc+nb);
    d_tau = d_bar_x(n+1);

    z = z + d_z;
    y = y + d_y;
    tau = tau + d_tau;

    Qz = Q*z;
    ATy = A'*y;
    Az = A*z;
    zTc = z'*c;
    yTb = y'*b;
    zTQz = z'*Qz;

    v = Qz - ATy +tau*c + gamma*r_z;
    w = Az - tau*b + gamma*r_y;
    kappa = -zTQz/tau - zTc + yTb + gamma*r_tau;  

    r_z = gamma*r_z;
    r_y = gamma*r_y;
    r_tau = gamma*r_tau;

    duality = z'*v+y'*w+tau*kappa;
    Dualitys = [Dualitys; duality];
    infeasibility = norm([r_z;r_y;r_tau]);
    Infeasibilitys = [Infeasibilitys; infeasibility];
    if (duality<=epsilon && infeasibility<=epsilon)
        break
    end
end
if(kappa<tau)
    status = "Optimal";
    z = z/tau;
else
    status = "Infeasible";
    z = z/tau;
end

