function MPO = BLBQMPO(J, Nsite,varargin)
    % Hubbard model with SU(2) chage, spin symmetry
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
    % near0 = 1e-40;
    
    [S, info] = getLocalSpace('Spin', 1/2);
    Ic = info.E;

    C = getIdentity(S,3,S,3);
    S2 = contract(S,2,S,1);
    SS = contract(S2,'24',C,'12');

    % Cd = permute(C,[1,2,3],'conj');
    % S2d = contract(Hconj(S),2,Hconj(S),1);
    % SSd = contract(S2d,'24',Cd,'12');

    SSd = Hconj(contract(S2,'24',C,'21'));

    Hloc = cell(4,4);
    Hloc{4,2} = S;
    Hloc{4,3} = J*SS;
    Hloc{4,4} = Ic;
    Hloc{1,1} = Ic;
    Hloc{2,1} = Hconj(S);
    Hloc{3,1} = SSd;

    % Hloc = cell(4,4);
    % Hloc{4,2} = I*S + near0*SS;
    % Hloc{4,3} = J*SS + near0*S;
    % Hloc{4,4} = Ic;
    % Hloc{1,1} = Ic;
    % Hloc{2,1} = Hconj(S) + near0*SSd;
    % Hloc{3,1} = SSd + near0*Hconj(S);

    MPOloc = getMPO(Hloc,Ic);

    MPO = QSpace(Nsite,1);
    MPO(1) = getMPO(Hloc,Ic, "start");
    MPO(end) = getMPO(Hloc,Ic, "end");

    % MPOloc = getMPOMK(Hloc,Ic);

    % MPO = QSpace(Nsite,1);
    % MPO(1) = getMPOMK(Hloc,Ic, "start");
    % MPO(end) = getMPOMK(Hloc,Ic, "end");
    
    for itN = 2:Nsite-1
        MPO(itN) = MPOloc;
    end
    
    
    if itagson 
        MPO = attachitags2MPO(MPO,'s','O');
    end    
    
    end

