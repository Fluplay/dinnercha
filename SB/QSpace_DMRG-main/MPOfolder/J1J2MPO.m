function MPO = J1J2MPO(J1, J2,Nsite,varargin)
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
    I = info.E;
    
    Hloc = cell(4,4);
    Hloc{1,1} = I;
    Hloc{2,1} = Hconj(S);
    Hloc{3,1} = 1e-40*Hconj(S);
    Hloc{end,1} = 1e-40*I;
    Hloc{end,end} = I;
    Hloc{3,2} = I;
    Hloc{end,2} = J1*S;
    Hloc{end,3} = J2*S;
    MPOloc = getMPO(Hloc,I);
    
    MPO = QSpace(Nsite,1);
    MPO(1) = getMPO(Hloc,I, "start");
    MPO(end) = getMPO(Hloc,I, "end");
    
    for itN = 2:Nsite-1
        MPO(itN) = MPOloc;
    end
    
    if itagson 
        MPO = attachitags2MPO(MPO,'s','O');
    end    
    
    end