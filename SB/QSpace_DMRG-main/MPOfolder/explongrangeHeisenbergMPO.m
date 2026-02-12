function MPO = explongrangeHeisenbergMPO(J, lambda, Nsite)

%% exponentially decaying long-range Heisenberg interaction

[S, info] = getLocalSpace('Spin', 1/2); % 1/2
I = info.E;

Hloc = cell(3,3);
Hloc{1,1} = I;
Hloc{2,1} = Hconj(S);
Hloc{end,2} = J*S;
Hloc{end,1} = 1e-40*I;
Hloc{end,end} = I;
Hloc{2,2} = lambda*I;

MPOloc = getMPOMK(Hloc,I);

MPO = QSpace(Nsite,1);
MPO(1) = getMPOMK(Hloc,I, "start");
MPO(end) = getMPOMK(Hloc,I, "end");

for itN = 2:Nsite-1
    MPO(itN) = MPOloc;
end

MPO = attachitags2MPO(MPO,'s','O');


end