function MPO = XXZMPO(delta,h, Nsite)

near0=1e-40;

[S, info] = getLocalSpace('Spin', 1/2, '-A');
I = info.E;

Sz = S(1);
Sm = S(2);
Sp = S(3);



Hloc = cell(5,5);
Hloc{1,1} = I;
Hloc{2,1} = Hconj(Sp);
Hloc{3,1} = Hconj(Sm);
Hloc{4,1} = Hconj(Sz);
Hloc{5,1} = h*Sz;
Hloc{5,2} = Sp;
Hloc{5,3} = Sm;
Hloc{5,4} = delta*Sz;
Hloc{5,5} = I;
MPOloc = getMPOMK(Hloc, I);

MPO = QSpace(Nsite,1);
MPO(1) = getMPOMK(Hloc,I, "start");
MPO(end) = getMPOMK(Hloc,I, "end");

for itN = 2:Nsite-1
    MPO(itN) = MPOloc;
end

MPO = attachitags2MPO(MPO,'s','O');



end