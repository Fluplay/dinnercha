function MPO = HubbardMPO(U, mu, t, Nsite)
    % Hubbard model with SU(2) chage, spin symmetry

    if U == 0
        error('please use realmin instead of 0')
    end
    if U == 0
        error('please use realmin instead of 0')
    end
    
    [F,Z,S,I] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1); %particle-hole symmetry
    Ic = I.E;
    nloc = contract(F,'13*',F,'13');


    Hloc = cell(4,4);
    Hloc{4,1} = U*contract((nloc-Ic),2,nloc,1)/2 - mu*nloc;
    Hloc{4,2} = t*F;
    Hloc{4,3} = -t*legflip(Hconj(F),3);
    Hloc{4,4} = Ic;
    Hloc{3,1} = legflip(contract(Z,2,F,1),3);
    Hloc{2,1} = contract(Z,2,Hconj(F),1);
    Hloc{1,1} = Ic;

    MPOloc = getMPOSB(Hloc,Ic);

    
    MPO = QSpace(Nsite,1);
    MPO(1) = getMPOSB(Hloc,Ic, "start");
    MPO(end) = getMPOSB(Hloc,Ic, "end");

    % MPOloc = getMPOMK(Hloc,Ic);

    
    % MPO = QSpace(Nsite,1);
    % MPO(1) = getMPOMK(Hloc,Ic, "start");
    % MPO(end) = getMPOMK(Hloc,Ic, "end");
    
    for itN = 2:Nsite-1
        MPO(itN) = MPOloc;
    end
    
    MPO = attachitags2MPO(MPO,'s','O');
    
    
    end
    