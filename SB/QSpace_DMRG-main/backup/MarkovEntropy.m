function [S, S1,S2] = MarkovEntropy(Md)

    Nsite = numel(Md);
    redMd = QSpace(Nsite,1);
    Dlr = QSpace(Nsite+2,1);

    S2 = zeros(1,Nsite); % first cell is empty
    S1 = zeros(1,Nsite); % first, last cell is empty

    for itN =1:Nsite
        redMd(itN) = tracePhy(Md(itN));
    end

    % Dlr(1) is empty
    % Dlr(end) is empty

    for itN = 1:Nsite-1
        if isempty(Dlr(itN))
            Dlr(itN+1) = redMd(itN);
        else
            Dlr(itN+1) = contract(Dlr(itN),1,redMd(itN),1);
        end
    end

    for itN = Nsite:-1:2
        if isempty(Dlr(itN-1))
            A = Md(itN-1);
        else
            A = contract(Dlr(itN-1),1,Md(itN-1),1);
        end
        if isempty(Dlr(itN+2))
            B = Md(itN);
        else
            B = contract(Md(itN),2,Dlr(itN+2),1);
        end
        % S2 section
        AB = contract(A,1,B,1,[1 3 2 4]);
        I1 = conj(getIdentity(AB,1,AB,2));
        AB = contract(I1,[1 2],AB,[1 2]);
        I2 = getIdentity(AB,2,AB,3);
        AB = contract(AB,[2 3],I2,[1 2]);
        S2(itN) = SEntropy(AB);
        % S1 section
        if itN ~= Nsite
            CB = contract(Dlr(itN),1,B,1);
            S1(itN) = SEntropy(CB);
        end
        S = sum(S2)-sum(S1);
        if isempty(Dlr(itN+2))
            Dlr(itN+1) = redMd(itN);
        else
            Dlr(itN+1) = contract(redMd(itN),Dlr(itN+2));
        end
    end
end

function T = tracePhy(T)
    if rank(T) == 4
        I = getIdentity(T,4);
        T = contract(T,[4 3],I,[1 2]);
    else
        I = getIdentity(T,3);
        T = contract(T,[3 2],I,[1 2]);
    end

end