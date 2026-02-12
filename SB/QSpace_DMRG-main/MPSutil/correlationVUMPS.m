function OO = correlationVUMPS(Al,Ar,C,O,n,rArray)

    % OO(i) = <O(n)O(n+rArray(i))>

    Nunit = numel(Al);
    rMax = max(rArray)+n-1;
    cellNum = ceil((rMax+1)/(Nunit*2));
    Mu = [repmat(Al,[2*cellNum 1]);contract(C(Nunit),2,Ar(1),1);Ar(2:Nunit)];%;repmat(Ar,[(cellNum -1) 1])];
    OO = zeros(1,length(rArray));
    for i = 1:length(rArray)
        r = rArray(i);
        OO(i) = correlationM12MPS(Mu,Mu,O,r,'Selective',n);
    end
end