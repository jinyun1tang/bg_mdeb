function dydt=sDEB(t,y)
%
%mDEB model for diatom growing on CO2 and inorganic N
%Jinyun Tang(jinyuntang@lbl.gov)
%History: April 12th, 2024.
%t: day,
%y: state variable
%
global par;
global funcg;
dydt=zeros(size(y));vid=par.vid;
N_conc=y(vid.N);  %inorganic N concentration, muM
C_conc=y(vid.C);  %dissolved CO2 concentration, muM
[I,h]=irad(t);        %irradiance muEm^-2 s^-1
xC=y(vid.xC);     %C reserve density, dimensionless
xN=y(vid.xN);     %N reserve density, dimensionless
mV=y(vid.mV);     %structural biomass concentration, muM

%reserve mobilization
jE0_C=par.vE0*xC;
jE0_N=par.vE0*xN;
%reserve assimilation
%N assimilation, mol N /mol mV /day
jA_N=par.jm_AN*N_conc/(N_conc+par.KN);
%photosynthesis, mol C/mol mV/day
j_ph=par.rho_PSU*I/(1./par.alpha+1./par.gamma*(1+par.beta/par.alpha)*I+...
    (par.beta/(par.gamma*par.delta)*I*I))*par.secs1day;

j_C=par.jm_AC*C_conc/(par.KC+C_conc);

jA_C=1./(0/par.jm_AC+1./j_C+1./j_ph-1./(j_ph+j_C));

%potential growth flux
jG_C=(jE0_C-par.jM_C)/par.yVC;
jG_N=(jE0_N-par.jM_N)/par.yVN;

%record respiration flux
%jX_C=min([par.jM_C,jE_C]);jX_N=min([par.jM_N,jE_N]);
%growth, rejection, respiration and senescence/lysis
if(jG_C > 0. && jG_N > 0.)
    %positive growth, solve using bisection iteration
    mu1=0;  
    mu2=1./(1./jG_C+1./jG_N-1./(jG_C+jG_N));
    f1=mu1-mu2;                      %function value at minimum growth rate < 0
    jG_C=(jE0_C-mu2*xC-par.jM_C)/par.yVC;
    jG_N=(jE0_N-mu2*xN-par.jM_N)/par.yVN;
    f2=mu2-SUK(jG_C,jG_N);           %function value at maximum growth rate > 0

    while(1)
        mu=(mu1+mu2)*0.5;        
        jG_C=(jE0_C-mu*xC-par.jM_C)/par.yVC;
        jG_N=(jE0_N-mu*xN-par.jM_N)/par.yVN;
        fnew=mu-SUK(jG_C,jG_N);
        if(abs(fnew)<1.e-6)
            break;
        end
        if(fnew>0)
            f2=fnew;
            mu2=mu;
        else
            f1=fnew;
            mu1=mu;
        end
    end


    jG_C=(jE0_C-mu*xC-par.jM_C)/par.yVC;
    jG_N=(jE0_N-mu*xN-par.jM_N)/par.yVN;
    
    jR_C=(jG_C-mu)*par.yVC;
    jR_N=(jG_N-mu)*par.yVN;
    jX_C=(par.yVC-par.nV_C)*mu+par.jM_C;
    jX_N=(par.yVN-par.nV_N)*mu+par.jM_N;
else
    
    %negative growth, structural biomass is lost/respired
    %to make up the maintenance cost
    mu=-max([-jG_C,-jG_N,0.]);
    jR_C=max([jG_C,0.])*par.yVC;
    jR_N=max([jG_N,0.])*par.yVN;
    jX_C=min([jE0_C,par.jM_C])-mu*par.nV_C;
    jX_N=min([jE0_N,par.jM_N])-mu*par.nV_N;
end

%update state variable trend
dydt(vid.xC)=jA_C-jE0_C+par.kapaEC*jR_C;
dydt(vid.xN)=jA_N-jE0_N+par.kapaEN*jR_N;
dydt(vid.mV)=(mu-h)*mV;
dydt(vid.C)=h*(par.Cr-C_conc)-(jA_C-jX_C)*mV;
dydt(vid.N)=h*(par.Nr-N_conc)-(jA_N-jX_N)*mV;

%dydt(vid.cbal)=(dydt(vid.xC)+mu*y(vid.xC)+mu*par.nV_C)*y(vid.mV)+dydt(vid.C)-par.h*(par.Cr-C_conc);
%dydt(vid.nbal)=(dydt(vid.xN)+mu*y(vid.xN)+mu*par.nV_N)*y(vid.mV)+dydt(vid.N)-par.h*(par.Nr-N_conc);
end


function [I,h]=irad(t)
global par;
h=par.h;
if t<=4
    I=250.;
    h=0.;
elseif t<=24.5
    I=180.*exp(3.*(cos(2.0.*pi.*(t+0.4))-1));
else
    I=100;
end

end

function mu=SUK(jGC,jGN)

mu=1./(1./jGC+1./jGN-1./(jGC+jGN));
end
