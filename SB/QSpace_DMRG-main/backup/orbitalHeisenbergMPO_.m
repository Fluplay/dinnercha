function MPO = orbitalHeisenbergMPO_(J, K, I, Nsite)
    % Hubbard model with SU(2) chage, spin symmetry
    near0 = 1e-40;
    
    [F,Z,S,IS]=getLocalSpace('FermionS','Acharge,SU2spin,SU2channel','NC',2);
    Ic = getsub(IS.E,[-1 1 1],2);
    % nloc = contract(F,'13*',F,'13');
    U = getIdentity(F,3,'-0');
    Eop = contract(getIdentity(U,1,U,2),2,U,'2*',[1 3 2]);
    x = getsub(Eop,[0 2 0],3);
    Sop = contract({F,'1*',F,1},'24',x,'21');
    x = getsub(Eop,[0 0 2],3); % select pseudo-spin operator
    Top = contract({F,'1*',F,1},'24',x,'21');
    x = getsub(Eop,[0 2 2],3); % select spin-orbit operator
    STop = contract({F,'1*',F,1},'24',x,'21') / 2;

    Sop = getsub(Sop,[-1 1 1], 2);
    Top = getsub(Top,[-1 1 1], 2);
    STop = getsub(STop,[-1 1 1], 2);

    % make dummy sectors
    % Sop_ = cat(Sop,Top*near0,STop*near0,3);
    % Top_ = cat(Sop*near0,Top,STop*near0,3);
    % STop_ = cat(Sop*near0,Top*near0,STop,3);
    % Hloc = cell(5,5);
    % Hloc{5,1} = near0*Ic;
    % Hloc{5,2} = J*Sop_;
    % Hloc{5,3} = K*Top_;
    % Hloc{5,4} = I*STop_;
    % Hloc{5,5} = Ic;
    % Hloc{1,1} = Ic;
    % Hloc{2,1} = Hconj(Sop_);
    % Hloc{3,1} = Hconj(Top_);
    % Hloc{4,1} = Hconj(STop_);

    Hloc = cell(5,5);
    % Hloc{5,1} = 1e5*contract(nloc-Ic,2,nloc-Ic,1);
    Hloc{5,2} = J*Sop;
    Hloc{5,3} = K*Top;
    Hloc{5,4} = I*STop;
    Hloc{5,5} = Ic;
    Hloc{1,1} = Ic;
    Hloc{2,1} = Hconj(Sop);
    Hloc{3,1} = Hconj(Top);
    Hloc{4,1} = Hconj(STop);

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

