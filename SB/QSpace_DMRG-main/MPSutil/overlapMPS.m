function a =  overlapMPS(M1,M2,Nsite)
    % <Describe>

    % Calculate overlap of M1 and M2

    % <Input>

    % M1, M2 : MPS
    % Nsite

    % <Output>

    % a : [QSpace tensor] 2-rank tensor. Trace of a is overlap.

    % Written by Subin Kim (01.02.2025)
    for itN = 1:Nsite

        if itN == 1
            Mval = updateLeft([],[],M1(itN),[],[],M2(itN));
        elseif itN == Nsite
            Mval = updateLeft(Mval,2,M1(itN),[],[],M2(itN));
        else
            Mval = updateLeft(Mval,2,M1(itN),[],[],M2(itN));
        end
    end
    % Mval = contract(Mval,2,M2(Nsite),1);

    % Mval = contract(Mval,[1 2 3],conj(M1(Nsite)), [1 2 3]);
    % trace(Mval)
    a = Mval;
end