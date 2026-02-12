function [EE, Sd, M, Mc, Mc0, Sd0,dwSet,kwSet] = entropyDividedMPS(M,Ds)
    Nsite = numel(M);
    Nc = round(Nsite/2);
    dwSet = zeros(1,Nsite);
    kwSet = zeros(1,Nsite);
    for itN = 1:Nc-1
        itag = getitags(M(itN), 2);
        if ~isempty(itag) && itag(end) == '*'
            itag = itag(1:end-1);
        end
        dimM = dim(M(itN));
        [M(itN),Sd,Vd] = svdSB(M(itN),[1 3],dimM(2),0,itag);
        M(itN) = permute(M(itN),[1 3 2]);
        V = contract(Sd,2,Vd,1);
        M(itN+1) = contract(V,2,M(itN+1),1,[1 3 4 2 5]);
        M(itN+1) = legfuse(M(itN+1),[4 5],'in');
        [M(itN+1),Sd,~,dw] = svdSB(M(itN+1),[1 2 3],Ds,0);
        dwSet(itN+1) = dw;
        kwSet(itN+1) = sqrt(trace(contract(Sd,2,Sd,'2*')));
        M(itN+1) = contract(M(itN+1),4,Sd,1);
    end
    for itN = Nsite:-1:Nc+2
        itag = getitags(M(itN), 1);
        if ~isempty(itag) && itag(end) == '*'
            itag = itag(1:end-1);
        end
        dimM = dim(M(itN));
        [U,Sd,M(itN)] = svdSB(M(itN),[1 4],dimM(1),0,itag);
        U = contract(U,3,Sd,1,[1 3 2]);
        M(itN-1) = contract(M(itN-1),2,U,1,[1 4 2 3 5]);
        M(itN-1) = legfuse(M(itN-1),[4 5],'in');
        [M(itN-1),Sd,~,dw] = svdSB(M(itN-1),[1 2 3],Ds,0);
        dwSet(itN-1) = dw;
        kwSet(itN-1) = sqrt(trace(contract(Sd,2,Sd,'2*')));
        M(itN-1) = contract(M(itN-1),4,Sd,1);
    end
    dimM = dim(M(Nc));
    [M(Nc),Sd,Vd] = svdSB(M(Nc),[1 3],dimM(2),0);
    Ml = contract(Sd,2,Vd,1);
    dimM = dim(M(Nc+1));
    [U,Sd,M(Nc+1)] = svdSB(M(Nc+1),[1 4],dimM(1),0);
    Mr = contract(U,3,Sd,1);
    Mc = contract(Ml,2,Mr,1,[1 4 2 3]);
    [Mc0, Sd0, ~] = svdSB(Mc,[1 2],1,0);
    [Mc, Sd, ~] = svdSB(Mc,[1 2],[],0);
    S2 = contract(Sd,2,Sd,'2*');
    % trace(S2)
    S2 = S2/trace(S2);
    EE = SEntropy(S2);
end