function MPO = orbitalHeisenbergAbelianMPO(J, K, I, f, Nsite)
    near0 = 1e-40;
    
    [~,~,S,Is]=getLocalSpace('FermionS','Acharge,Aspin','NC',1);
    Ic = Is.E;
    Ic = getsub(Ic,[0 1],2) + getsub(Ic,[0 -1],2);

    Ic1 = addSymmetry(Ic,'A','q',1);
    Ic2 = addSymmetry(Ic,'A','q',0);
    Ic3 = addSymmetry(Ic,'A','q',-1);
    Ic = Ic1 + Ic2 + Ic3;

    Sop = QSpace(1,3);
    Sop(1) = addSymmetry(S(1),'A','q',1) + addSymmetry(S(1),'A','q',-1);
    Sop(2) = addSymmetry(S(2),'A','q',1) + addSymmetry(S(2),'A','q',-1);
    Sop(3) = addSymmetry(S(3),'A','q',1) + addSymmetry(S(3),'A','q',-1);
    Sz = Sop(3);
    Sop(3) = appendSingletons(Sop(3),'  -');
    Top = QSpace(1,3);
    Top(1) = addSymmetry(S(1),'A','pos',2,'q',1) + addSymmetry(S(1),'A','pos',2,'q',-1);
    Top(2) = addSymmetry(S(2),'A','pos',2,'q',1) + addSymmetry(S(2),'A','pos',2,'q',-1);
    Top(3) = addSymmetry(S(3),'A','pos',2,'q',1) + addSymmetry(S(3),'A','pos',2,'q',-1);
    Tz = Top(3);
    Top(3) = appendSingletons(Top(3),'  -');

    Stot = Sop(1) + Sop(2) + Sop(3);
    Ttot = Top(1) + Top(2) + Top(3);
    STtot = contract(Stot,2,Ttot,1,[1 3 2 4]);
    C = getIdentity(STtot,3,STtot,4);
    STtot = contract(STtot,[3 4],C,[1 2]);

    Stot = Stot + near0*Ttot + near0*STtot;
    Ttot = near0*Stot + Ttot + near0*STtot;
    STtot = near0*Stot + near0*Ttot + STtot;

    % make dummy sectors
    Hloc = cell(5,5);
    if f > 0
        Hloc{5,1} =  -f*Sz - f*Tz;
    end
    Hloc{5,2} = J*Stot;
    Hloc{5,3} = K*Ttot;
    Hloc{5,4} = I*STtot;
    Hloc{5,5} = Ic;
    Hloc{1,1} = Ic;
    Hloc{2,1} = Hconj(Stot);
    Hloc{3,1} = Hconj(Ttot);
    Hloc{4,1} = Hconj(STtot);


    % MPOloc = getMPOSB(Hloc,Ic);

    
    % MPO = QSpace(Nsite,1);
    % MPO(1) = getMPOSB(Hloc,Ic, "start");
    % MPO(end) = getMPOSB(Hloc,Ic, "end");

    % MPOloc = getMPOMK(Hloc,Ic);

    
    % MPO = QSpace(Nsite,1);
    % MPO(1) = getMPOMK(Hloc,Ic, "start");
    % MPO(end) = getMPOMK(Hloc,Ic, "end");
    
    MPOloc = getMPO(Hloc,Ic);

    
    MPO = QSpace(Nsite,1);
    MPO(1) = getMPO(Hloc,Ic, "start");
    MPO(end) = getMPO(Hloc,Ic, "end");
    
    for itN = 2:Nsite-1
        MPO(itN) = MPOloc;
    end
    
    MPO = attachitags2MPO(MPO,'s','O');
    
    
    end

