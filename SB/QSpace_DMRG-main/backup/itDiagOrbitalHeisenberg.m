function [Minit, EGSiter, Eeig] = itDiagOrbitalHeisenberg(J,K,I, Nsite, Nkeep)
    % MPS diagonalization not with MPO, but with neighboor interaction.
    [F,Z,S,IS]=getLocalSpace('FermionS','Acharge,SU2spin,SU2channel','NC',2);
    Ic = getsub(IS.E,[-1 1 1],2);
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
    
    Minit = QSpace(Nsite,1);
    
    tobj = tic2;
    
    vac = getvac(Ic);
    Aprev = vac;
    Hprev = getIdentity(vac,2);
    % Aprev.info.itags = {'start','start*'};
    
    vacQ = vac.Q{1};
    
    for itN = 1:Nsite

        Anow = getIdentity(Aprev,2,Ic,2,[1 3 2]);
        if itN == 1
            Hnow = 0*updateLeft(Hprev, 2, Anow, [],[] , Anow);
        else
            Hnow = updateLeft(Hprev, 2, Anow, [],[] , Anow);
        end


        if itN > 1
            HS = updateLeft(Sprev, 3, Anow, Hconj(Sop), 3, Anow);
            HT = updateLeft(Tprev, 3, Anow, Hconj(Top), 3, Anow);
            HST = updateLeft(STprev, 3, Anow, Hconj(STop), 3, Anow);

            Hnow = Hnow + HS + HT + HST;
        end

        if itN ~= Nsite
            [Eeig, Ieig] = eigQS((Hnow+Hnow')/2, 'Nkeep', Nkeep);
            Ieig.AK = QSpace(Ieig.AK);
            Aprev = contract(Anow,2,Ieig.AK,1,[1 3 2]);
            Hprev = contract(conj(Ieig.AK), 1, contract(Hnow, 2, Ieig.AK, 1), 1);
            Minit(itN) = Aprev;

        else
            [Eeig, Ieig] = eigQS((Hnow+Hnow')/2, 'Nkeep', 1);
            Ieig.AK = QSpace(Ieig.AK);
            Aprev = contract(Anow,2,Ieig.AK,1,[1 3 2]);
            Hprev = contract(conj(Ieig.AK), 1, contract(Hnow, 2, Ieig.AK, 1), 1);
            Aprev = legflip(Aprev,2)
            Minit(itN) = Aprev;
        end


        Sprev = J*updateLeft([],[],Aprev,Sop,3,Aprev);
        Tprev = K*updateLeft([],[],Aprev,Top,3,Aprev);
        STprev = I*updateLeft([],[],Aprev,STop,3,Aprev);
        
    end
    EGSiter = Eeig(1);
    
    toc2(tobj,'-v');
    
    end