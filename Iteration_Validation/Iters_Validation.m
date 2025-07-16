clear,close,clc
rng(20250605)
nc = 20;
nb = 4*nc;
n = nc + nb;
eta = 0.414213/sqrt(n+1);
gamma = 1-eta;
theory_iters_1em6 = ceil(log((n+1)/1e-6)/(-log(gamma)));
theory_iters_1em8 = ceil(log((n+1)/1e-8)/(-log(gamma)));

%%
num_cond = 1e3;
v = [num_cond; (num_cond-1)*rand(nc-2,1)+1; 1];
[U,~] = qr(rand(nc,nc));
Q = U*diag(v)*U';
c = ones(nc,1);
A = full(2*sprand(nb,nc,0.15));
b = ones(nb,1);
[Iters,Dualitys,Infeasibilitys] = EIQP_Iters(Q,c,A,b,1e-8);

figure(1)
semilogy(Iters,Dualitys,'b--','LineWidth',2)
hold on 
semilogy(Iters,Infeasibilitys,'r--','LineWidth',2)
xlabel('Number of Iterations')
xline(theory_iters_1em6,'c--','LineWidth', 2)
xline(theory_iters_1em8,'g--','LineWidth', 2)
legend()

%%
nc = 40;
nb = 4*nc;
n = nc + nb;
eta = 0.414213/sqrt(n+1);
gamma = 1-eta;
theory_iters_1em6 = ceil(log((n+1)/1e-6)/(-log(gamma)));
theory_iters_1em8 = ceil(log((n+1)/1e-8)/(-log(gamma)));

%%
num_cond = 1e3;
v = [num_cond; (num_cond-1)*rand(nc-2,1)+1; 1];
[U,~] = qr(rand(nc,nc));
Q = U*diag(v)*U';
c = ones(nc,1);
A = full(2*sprand(nb,nc,0.15));
b = ones(nb,1);
[Iters,Dualitys,Infeasibilitys] = EIQP_Iters(Q,c,A,b,1e-8);


semilogy(Iters,Dualitys,'b','LineWidth',2)
semilogy(Iters,Infeasibilitys,'r','LineWidth',2)
xlabel('Number of Iterations')
xline(theory_iters_1em6,'c','LineWidth', 2)
xline(theory_iters_1em8,'g','LineWidth', 2)
yline(1e-6,'k--','LineWidth', 2)
yline(1e-8,'k--','LineWidth', 2)
legend('Duality gap (n=100)','Infeasibility residual norm (n=100)','Theoretical exact iterations (n=100,\epsilon=10^{-6}))','Theoretical exact iterations (n=100,\epsilon=10^{-8}))','Duality gap (n=200)','Infeasibility residual norm (n=200)','Theoretical exact iterations (n=200,\epsilon=10^{-6}))','Theoretical exact iterations (n=200,\epsilon=10^{-8}))')
