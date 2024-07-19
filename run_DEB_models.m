close all;
clear all;
%run and evaluate both DEB models
%addpath(genpath('/Users/jinyuntang/work/github/matlab_tools'))
%%
%configure models
global par;
par=setpar();

y0=zeros(1,par.nvars);
y0(par.vid.C)=1.e3;     %uM, Mx1.e-3, 1.e-6*1.e3mol/m3=1.e-3 mol/m3
y0(par.vid.N)=70;       %uM, 1.e-3 mol/m3
y0(par.vid.xC)=0.1;     %C-reserve density
y0(par.vid.xN)=0.01;    %N-reserve density
y0(par.vid.mV)=100.;    %uM, 1.e-3 mol/m3

tspan=[0,30];
ns=3200;
%%
%run models
tic;
[t1,y1]=euler_driver_sdeb(y0,tspan,ns);
elapsedTime1 = toc;

[t2,y2]=euler_driver_mdeb(y0,tspan,ns);
elapsedTime2 = toc;
%%
%print timing info
fprintf('sdeb time=%f seconds\n',elapsedTime1);
fprintf('mdeb time=%f seconds\n',elapsedTime2-elapsedTime1);
fprintf('sDEB is %f%% slower than mDEB\n',(elapsedTime1/(elapsedTime2-elapsedTime1)-1)*100);
%%
%compute the radiation forcing
t=(0:0.1:30);

I=(t<=4).*250+(t>4 & t<=24.5).*180.*exp(3.*(cos(2.0.*pi.*(t+0.4))-1))...
    +(t>24.5).*100;

%%
%plot the results
load('C_conc.mat');
load('N_conc.mat');

fig=1;
figure(fig);
axs=zeros(3,1);
axs(1)=subplot(3,1,1);
plot(t,I,'LineWidth',2);xlim([20,30]);
ylabel('PFD ($\mu$E m$^{-2}$ s$^{-1}$)', 'Interpreter', 'latex');

axs(2)=subplot(3,1,2);
plot(t1,(y1(:,par.vid.xC)+1).*y1(:,par.vid.mV).*1.e-3,'b','LineWidth',2);
hold on;
plot(t2,(y2(:,par.vid.xC)+1).*y2(:,par.vid.mV).*1.e-3,'k','LineWidth',2);

plot(C_conc(:,1),C_conc(:,2),'r+','MarkerSize',12);xlim([20,30]);
ylabel('POC concentration (M) $\times 10^{-3}$', 'Interpreter', 'latex');

axs(3)=subplot(3,1,3);
h(1)=plot(t1,(y1(:,par.vid.xN)+par.nV_N).*y1(:,par.vid.mV).*1.e-2,'b','LineWidth',2);hold on;
h(2)=plot(t2,(y2(:,par.vid.xN)+par.nV_N).*y2(:,par.vid.mV).*1.e-2,'k','LineWidth',2);
h(3)=plot(N_conc(:,1),N_conc(:,2),'r+','MarkerSize',12);xlim([20,30]);
legend(h,'sDEB model','mDEB model','Measurement');
ylabel('PON concentration (M) $\times 10^{-4}$', 'Interpreter', 'latex');
xlabel('Day');
set(axs,'FontSize',16);
set(axs(1:2), 'XTickLabel', []);
set(axs(3), 'Position', [0.12, 0.08, 0.8, 0.255]);
set(axs(2), 'Position', [0.12, 0.4, 0.8, 0.255]);
set(axs(1), 'Position', [0.12, 0.72, 0.8, 0.255]);
set(fig,'color','w');
%%
%put figure tag
%can be commented out
put_tag(fig,axs(1),[.025,.9],'(a)',16);
put_tag(fig,axs(2),[.025,.9],'(b)',16);
put_tag(fig,axs(3),[.025,.9],'(c)',16);
