function MPO = BLBQAbelianMPO(J, I, Nsite)
    % Hubbard model with SU(2) chage, spin symmetry
    near0 = 1e-40;
    
    [S, info] = getLocalSpace('Spin',1,'-A');
    Ic = info.E;

    Sz = S(1);
    vac = getvac(Sz,'-1d');
    Sz = contract(Sz,3,vac,1);
    Hloc = cell(14,14);
    Hloc{14,1} = 1e-6*Sz; %for breaking symmetry
    Hloc{14,14} = Ic;
    Hloc{1,1} = Ic;
    for i = 2:4
        Hloc{14,i} = I*S(i-1);
        Hloc{i,1} = Hconj(S(i-1));
    end
    for i = 0:8
        Hloc{14,i+5} = J*makeSS(S(floor(i/3)+1),S(mod(i,3)+1));
        % Hloc{i+5,1} = Hconj(makeSS(S(floor(i/3)+1),S(mod(i,3)+1)));
        Hloc{i+5,1} = makeSSd(S(floor(i/3)+1),S(mod(i,3)+1));
    end
    MPOloc = getMPOMK(Hloc,Ic);

    MPO = QSpace(Nsite,1);
    MPO(1) = getMPOMK(Hloc,Ic, "start");
    MPO(end) = getMPOMK(Hloc,Ic, "end");
    
    % MPOloc = getMPOSB(Hloc,Ic);

    % MPO = QSpace(Nsite,1);
    % MPO(1) = getMPOSB(Hloc,Ic, "start");
    % MPO(end) = getMPOSB(Hloc,Ic, "end");

    for itN = 2:Nsite-1
        MPO(itN) = MPOloc;
    end
    
    MPO = attachitags2MPO(MPO,'s','O');
    
    
    end
function SS = makeSS(S1,S2)


    C = getIdentity(S1,3,S2,3);
    S2 = contract(S1,2,S2,1);
    SS = contract(S2,'24',C,'12');
end
function SSd = makeSSd(S1,S2)
    C = getIdentity(S1,3,S2,3);
    Cd = permute(C,[1,2,3],'conj');
    S2d = contract(Hconj(S1),2,Hconj(S2),1);
    SSd = contract(S2d,'24',Cd,'12');
end
