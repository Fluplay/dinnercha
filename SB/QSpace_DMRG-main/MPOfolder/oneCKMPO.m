function MPO = oneCKMPO(t, J, TBchainlength)

if J == 0
    error('please use realmin instead of 0')
end
if t == 0
    error('please use realmin instead of 0')
end

[F, Z, Sc, Ic] = getLocalSpace('FermionS', 'Acharge, SU2spin');
Ic = Ic.E;
[Simp, Iimp] = getLocalSpace('Spin', 1/2);
Iimp = Iimp.E;
Simp = addSymmetry(Simp, 'A', 'pos' ,1);
Iimp = addSymmetry(Iimp, 'A', 'pos' ,1);


Hlocimp = cell(1,3);
Hlocimp{1,1} = realmin*Iimp;
Hlocimp{1,2} = J*Simp;
Hlocimp{1,3} = Iimp;
MPOimp = getMPOMK(Hlocimp,Iimp, 'start');

Hloccimp = cell(3,4);
Hloccimp{1,1} = Ic;
Hloccimp{2,1} = Hconj(Sc);
Hloccimp{3,1} = Ic*realmin;
Hloccimp{3,2} = t*legflip(Hconj(F),3);
Hloccimp{3,3} = t*F;
Hloccimp{3,4} = Ic;
MPOcimp = getMPOMK(Hloccimp, Ic);

Hlocc = cell(4,4);
Hlocc{1,1} = Ic;
Hlocc{2,1}  =legflip(contract(Z,2,F,1),3);
Hlocc{3,1} = contract(Hconj(F),2, Z', 1, [1 3 2]);
Hlocc{4,1} = Ic*realmin;
Hlocc{4,2} = t*legflip(Hconj(F),3);
Hlocc{4,3} = t*F;
Hlocc{4,4} = Ic;
MPOc = getMPOMK(Hlocc, Ic);


MPO = QSpace(TBchainlength+1,1);
MPO(1) = MPOimp;
MPO(2) = MPOcimp;
MPO(end) = getMPOMK(Hlocc, Ic, 'end');

for itN = 3:TBchainlength
    MPO(itN) = MPOc;
end

MPO = attachitags2MPO(MPO,'s','O');

end
