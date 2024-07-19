function par=setpar()
%set parameter for the mDEB model
%reference:
%Lorena et al. (2010), Stylized facts in microalgal growth: interpretation in a dynamic energy budget context
par.vid.C=1;
par.vid.N=2;
par.vid.xC=3;
par.vid.xN=4;
par.vid.mV=5;
par.nvars=par.vid.mV;
%par.vid.cbal=7;
%par.vid.nbal=8;
par.vE=0.55;           %d^-1, reserve turnover rate for mDEB, tunable
par.vE0=0.95;           %d^-1, reserve turnover rate for sDEB
par.KC=3.20;           %uM, half-saturation concentration for C, K=v/D, D=v/K
par.KN=0.43;           %uM, half-saturation concentration for N
par.alpha=0.0019;      %1/(mmolPSU muEm^-2), PSU excitement coefficient
par.beta=5.8e-7;       %1/(mmolPSU muEm^-2), PSU inhibition coefficient
par.gamma=0.1460;      %mmolPSU^-1 s^-1, PSU relaxation rate
par.delta=4.8e-4;      %mmolPSU mol mV^-1, PSU recovery rate
par.rho_PSU=0.365e-3;  %mmolPSU mol mV^-1, PSU density per mol structural biomass
par.secs1day=86400;    %seconds per day
par.jm_AN=1.0e-0;      %molN mol mV^-1 d^-1, maximum volume-specific N assimilation
par.jm_AC=5.1e-0;      %molC mol mV^-1 d^-1, maximum volume-specific C assimilation
par.jM_C=0.054;        %molEC mol mV^-1 d^-1, volume-specific maintenance cost paid by C-reserve
par.jM_N=0.012;        %molEN mol mV^-1 d^-1, volume-specific maintenance cost paid by N-reverve
par.nV_C=1;            %molC /mol mV, mol C per mol of structural C biomass
par.nV_N=0.04;         %molN /mol mV, mol N per mol of structural C biomass
par.yVC=1.25;          %mol mV/mol C, reserve carbon amount required for producing one unit of structural biomass
par.yVN=0.04;          %mol mV/mol N, N use efficiency for converting N reserve to structural biomass
par.h=0.4;             %1/day, dilution rate
par.kapaEC=0.7;        %fraction of rejected reserve reincoporateedd into reserve, aka remobilization efficiency
par.kapaEN=0.7;        %fraction of rejected reserve reincoporateedd into reserve, aka remobilization efficiency
par.Cr=1.e3;           %supply concentration of CO2, uM
par.Nr=70;             %supply concentration inorganic N, uM
end
