function MPO = orbitalHeisenbergPeriodicMPO(J, K, I, Nsite)
    
    [F,~,~,IS]=getLocalSpace('FermionS','Acharge,SU2spin,SU2channel','NC',2);
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
    MPO = QSpace(Nsite,1);

    % Make Ic 4rank
    IcS = getIdentity(Sop,3);
    IcS = directprod(Ic,IcS);
    IcT = getIdentity(Top,3);
    IcT = directprod(Ic,IcT);
    IcST = getIdentity(STop,3);
    IcST = directprod(Ic,IcST);


    for itN = 1:Nsite
        Hloc = cell(8,8);
        Hloc{8,2} = J*Sop;
        Hloc{8,3} = K*Top;
        Hloc{8,4} = I*STop;
        Hloc{8,8} = Ic;
        Hloc{1,1} = Ic;
        Hloc{2,1} = Hconj(Sop);
        Hloc{3,1} = Hconj(Top);
        Hloc{4,1} = Hconj(STop);
        if itN == 1
            Hloc{8,5} = J*Sop;
            Hloc{8,6} = K*Top;
            Hloc{8,7} = I*STop;
        elseif itN == Nsite
            Hloc{5,1} = Hconj(Sop);
            Hloc{6,1} = Hconj(Top);
            Hloc{7,1} = Hconj(STop);
        else
            Hloc{5,5} = IcS;
            Hloc{6,6} = IcT;
            Hloc{7,7} = IcST;
        end

        if itN == 1 
            MPO(1) = getMPO(Hloc,Ic, "start");
        elseif itN == Nsite
            MPO(end) = getMPO(Hloc,Ic, "end");
        else
            MPO(itN) = getMPO(Hloc,Ic);
        end
    end
    
    
    MPO = attachitags2MPO(MPO,'s','O');
    
    
    end