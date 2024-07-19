function plot_obs()
load('C_conc.mat');
load('N_conc.mat');

t=(0:0.05:30);
par=setpar();
I=(t<=4).*250+(t>4 & t<=24).*180.*exp(3.*(cos(2.0.*pi.*(t+0.4))-1))...
    +(t>24).*100;

subplot(2,1,1);
plot(C_conc(:,1),C_conc(:,2),'+');
subplot(2,1,2);
plot(N_conc(:,1),N_conc(:,2),'+');

j_ph=par.rho_PSU.*I./(1./par.alpha+1./par.gamma*(1+par.beta/par.alpha).*I+...
    (par.beta/(par.gamma*par.delta).*I.*I))*par.secs1day; %1/day

figure;
plot(t,j_ph);
xlim([20,30]);
end




