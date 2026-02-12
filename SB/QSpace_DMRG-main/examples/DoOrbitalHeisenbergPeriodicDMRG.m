%% Operators
[S, info] = getLocalSpace('Spin', 1/2);
Sop = addSymmetry(S,'SU2','q',1);
Top = addSymmetry(S,'SU2','q',1,'pos',1);
STop = contract(Sop,2,Top,1,[1 3 2 4]);
C = getIdentity(STop,3,STop,4);
STop = contract(STop,[3 4],C,[1 2]);
Ic = getIdentity(STop,1);
%% Settings
J = 0.7; %0.25
K = 1; % -0.2
I = 1;
f = 0;
Nsite = 40;
MPO = orbitalHeisenbergPeriodicNoChargeMPO(J, K, I, Nsite);
%% MPS initialize
Nkeepinit=200;
[Minit, Egs, Eeig] = MPSinitialize(MPO, Nkeepinit,'E','Enum',1);
%% DMRG 2site-N40
Minit = MPSinitializeQn(MPO,0,0,'noCharge');
Nsweep = 10;
NkeepDMRG = 80;
alpha = 1.05;
[MDMRG2, Egs,Eiter,Sv] = DMRG_GS_2site_alpha_QSpace (Minit,MPO,NkeepDMRG,Nsweep,alpha,'raiseFlat',[1,0],'maxNkeep',1000);
EEspec(Sv)
correlationMPS(MDMRG2,Sop,1)
%% Spin-Orbital Entanglement
% reduce bond dimension for spin-orbital entanblement analysis
NkeepReduce = 70;
NsweepReduce = 1;
delta = 0.2;
alpha = 1.05;
Mr = DMRG_GS_RSVD_QSpace (MDMRG2,MPO,NkeepReduce,NsweepReduce,delta,alpha,'Canonical',true);
% center canonical form
[Mr,Sc] = canonFormMK(Mr,round(Nsite/2),[]);
M = divideOrbitalLeg(Mr,'noCharge');
Ds=90;
[EE,Sd,Mspin, Mcm, Mc0, Sd0, dwSet,kwSet,traceS2] = entropyDividedCanMPS(M,Sc,Ds);