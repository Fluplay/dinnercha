function [EE, Sd, Mc, Mc0, Sd0,dwSet,kwSet,traceS2,SdAl,SdbefAl,SdAr,SdbefAr] = entropyuMPS(Al,Ar,C,Ds, tol)

    % C is a single tensor
    % Written by Subin Kim (01.02.2025)
    iterMax = 50;
    Nunit = numel(Al);
    dwSet = zeros(1,iterMax);
    kwSet = zeros(1,iterMax);
    SdbefAl = QSpace(Nunit,1);
    SdAl = QSpace(Nunit,1);
    itN = 0;
    AlS = Al(1);
    converge = false;

    itN = itN+1;
    dimM = dim(AlS);
    [~,Sd,Vd] = svdSB(AlS,[1 3],dimM(2),0);
    V = contract(Sd,2,Vd,1);
    AlS = contract(V,2,Al(circularIndex(itN+1,Nunit)),1,[1 3 4 2 5]);
    AlS = legfuse(AlS,[4 5],'in');
    [AlS,Sd,~,dw] = svdSB(AlS,[1 2 3],Ds,0);
    dwSet(itN) = dw;
    kwSet(itN) = sqrt(trace(contract(Sd,2,Sd,'2*')));
    AlS = contract(AlS,4,Sd,1);
    SdbefAl(circularIndex(itN+1,Nunit)) = SdAl(circularIndex(itN+1,Nunit));
    SdAl(circularIndex(itN+1,Nunit)) = Sd;
    while ~converge && iterMax > itN
        for rep = 1:2
            itN = itN+1;
            dimM = dim(AlS);
            [~,Sd,Vd] = svdSB(AlS,[1 3],dimM(2),0);
            V = contract(Sd,2,Vd,1);
            AlS = contract(V,2,Al(circularIndex(itN+1,Nunit)),1,[1 3 4 2 5]);
            AlS = legfuse(AlS,[4 5],'in');
            [AlS,Sd,~,dw] = svdSB(AlS,[1 2 3],Ds,0);
            dwSet(itN) = dw;
            kwSet(itN) = sqrt(trace(contract(Sd,2,Sd,'2*')));
            AlS = contract(AlS,4,Sd,1);
            SdbefAl(circularIndex(itN+1,Nunit)) = SdAl(circularIndex(itN+1,Nunit));
            SdAl(circularIndex(itN+1,Nunit)) = Sd;
            try
                if(max(norm(SdbefAl-SdAl))) <  tol
                    converge = true;
                end
            catch
                SdbefAl
                SdAl
            end
        end
    end

    SdbefAr = QSpace(Nunit,1);
    SdAr = QSpace(Nunit,1);
    itN = 1;
    ArS = Ar(2);
    converge = false;

    itN = itN-1;
    dimM = dim(ArS);
    [U,Sd,~] = svdSB(ArS,[1 4],dimM(1),0);
    U = contract(U,3,Sd,1,[1 3 2]);
    ArS = contract(Ar(circularIndex(itN-1,Nunit)),2,U,1,[1 4 2 3 5]);
    ArS = legfuse(ArS,[4 5],'in');
    [ArS,Sd,~,dw] = svdSB(ArS,[1 2 3],Ds,0);
    dwSet(-itN+1) = dw;
    kwSet(-itN+1) = sqrt(trace(contract(Sd,2,Sd,'2*')));
    ArS = contract(ArS,4,Sd,1);
    SdbefAr(circularIndex(itN-1,Nunit)) = SdAr(circularIndex(itN-1,Nunit));
    SdAr(circularIndex(itN-1,Nunit)) = Sd;
    while ~converge  && iterMax > -itN+1
        for rep = 1:2
            itN = itN-1;
            dimM = dim(ArS);
            [U,Sd,~] = svdSB(ArS,[1 4],dimM(1),0);
            U = contract(U,3,Sd,1,[1 3 2]);
            ArS = contract(Ar(circularIndex(itN-1,Nunit)),2,U,1,[1 4 2 3 5]);
            ArS = legfuse(ArS,[4 5],'in');
            [ArS,Sd,~,dw] = svdSB(ArS,[1 2 3],Ds,0);
            dwSet(-itN+1) = dw;
            kwSet(-itN+1) = sqrt(trace(contract(Sd,2,Sd,'2*')));
            ArS = contract(ArS,4,Sd,1);
            SdbefAr(circularIndex(itN-1,Nunit)) = SdAr(circularIndex(itN-1,Nunit));
            SdAr(circularIndex(itN-1,Nunit)) = Sd;
            try
                if(max(norm(SdbefAr-SdAr))) <  tol
                    converge = true;
                end
            catch
                SdbefAr
                SdAr
            end
        end
    end
    dimM = dim(AlS);
    [~,Sd,Vd] = svdSB(AlS,[1 3],dimM(2),0);
    Ml = contract(Sd,2,Vd,1);
    dimM = dim(ArS);
    [U,Sd,~] = svdSB(ArS,[1 4],dimM(1),0);
    Mr = contract(U,3,Sd,1);
    Ml = contract(Ml,2,C,1,[1 3 2]);
    Mc = contract(Ml,2,Mr,1,[1 4 2 3]);
    [Mc0, Sd0, ~] = svdSB(Mc,[1 2],1,0);
    [Mc, Sd, ~] = svdSB(Mc,[1 2],[],0);
    S2 = contract(Sd,2,Sd,'2*');
    traceS2 = trace(S2);
    S2 = S2/trace(S2);
    EE = SEntropy(S2);
end
