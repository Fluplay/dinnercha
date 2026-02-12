function MPO = orbitalHeisenbergPeriodicNoChargeMPO(J, K, I, Nsite,varargin)
    itagson = true;
    while ~isempty(varargin)
        switch varargin{1}
            case 'itagsoff'
                itagson = false;
                varargin(1) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end
    [S, info] = getLocalSpace('Spin', 1/2);
    Sop = addSymmetry(S,'SU2','q',1);
    Top = addSymmetry(S,'SU2','q',1,'pos',1);
    STop = contract(Sop,2,Top,1,[1 3 2 4]);
    C = getIdentity(STop,3,STop,4);
    STop = contract(STop,[3 4],C,[1 2]);
    Ic = getIdentity(STop,1);
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
        Hloc{8,8} = Ic;
        Hloc{1,1} = Ic;
        Hloc{2,1} = Hconj(Sop);
        Hloc{3,1} = Hconj(Top);
        Hloc{4,1} = Hconj(STop);

        if itN == 1 
            MPO(1) = getMPO(Hloc,Ic, "start");
        elseif itN == Nsite
            MPO(end) = getMPO(Hloc,Ic, "end");
        else
            MPO(itN) = getMPO(Hloc,Ic);
        end
    end
    
    
    if itagson 
        MPO = attachitags2MPO(MPO,'s','O');
    end
    
    
    end