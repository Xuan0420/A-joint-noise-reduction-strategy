function SeMin=objfun(para,x)
% para(1)-alpha penalty parameter
% para(2)-K number of components
% x=Noise signal
alpha=round(para(1));
K=round(para(2));

tau = 0;            % noise-tolerance (no strict fidelity enforcement)
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-7;

%% VMD
[u, u_hat, omega] = VMD(x, alpha, tau, K, DC, init, tol);
%% Calculate envelope entropy
[m,n]=size(u);
for i=1:m
    En(i)=BaoluoEntropy(u(i,:));
end
%% Minimum local envelope entropy   
SeMin=min(En);
