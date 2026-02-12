function Hval = HvalMPS(M,MPO,Nsite)
    % <Description>

    % Calculate Hamiltonian expected value for MPS M.

    % <Input>

    % M : MPS
    % MPO : MPO
    % Nsite

    % <Output>

    % Hval : [QSpace tensor] 2-rank tensor. Contracted MPS and MPO. Trace of Hval is Hamiltonian expected value.

    % Written by Subin Kim (01.02.2025)

    %% Energy eigenvalue
    vac = getvac(MPO(1),'-1d');
    Hval = [];
    for itN = 1:Nsite
        MPOsite = MPO(itN);
        if itN == 1
            MPOsite = contract(legflip(vac,1),1,MPOsite,3);
        elseif itN == Nsite
            MPOsite = contract(MPOsite,4,vac, 1);
        end

        if itN == 1
            Hval = updateLeft([],[],M(itN),MPOsite,3,M(itN));
        elseif itN == Nsite
            Hval = updateLeft(Hval,3,M(itN),MPOsite,3,M(itN));
        else
            Hval = updateLeft(Hval,3,M(itN),MPOsite,4,M(itN));
        end
    end
end