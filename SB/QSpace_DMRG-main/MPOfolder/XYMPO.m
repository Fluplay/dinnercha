function MPO = XYMPO(J, Nsite)

[S, info] = getLocalSpace('Spin', 1/2, '-A');
I = info.E;

Sz = S(1);
Sm = S(2);
Sp = S(3);

Hloc = cell(4,4);
Hloc{1,1} = I;
Hloc{2,1} = Hconj(Sp);
Hloc{3,1} = Hconj(Sm);
Hloc{4,1} = realmin*I;
Hloc{4,2} = J*Sp;
Hloc{4,3} = J*Sm;
Hloc{4,4} = I;
MPOloc = getMPOMK(Hloc, I);

MPO = QSpace(Nsite,1);
MPO(1) = getMPOMK(Hloc,I, "start");
MPO(end) = getMPOMK(Hloc,I, "end");

for itN = 2:Nsite-1
    MPO(itN) = MPOloc;
end

MPO = attachitags2MPO(MPO,'s','O');



end