function [HLout, HLin, HRin, HRout, resNormHL, resNormHR] = HeffTerms(Al, Ar, C, h, itN,nKrylov,tol,iterNum)
    % < Description >

    % [HLout, HLin, HRin, HRout] = HeffTerms(Al, Ar, h, itN)


    % HLout : ...*TL(itN-1)
    % HLin : ...*TL(itN)
    % HRin : TR(itN)* ...
    % HRout : TR(itN+1)* ...
    % resNormHL(nunit,iter)
    % resNormHR(nunit,iter)
    

    Nunit = numel(Al);

    HLoutArr = QSpace(Nunit,1);
    HLinArr = HLoutArr;
    % Calculate 
    for n = 1:Nunit
        [HL,resNormHL(n,1)] =  gmres_HL_QSpace(Al,C,h,n,nKrylov,tol);
        % resNormHL(n,1)
        for iter = 1:iterNum-1
            [HL,resNormHL(n,iter+1)] = gmres_HL_QSpace(Al,C,h,n,nKrylov,tol,'x0',HL);
            % resNormHL(n,iter+1)
        end

        if circularIndex(n+1,Nunit) == itN
            HLinArr(n) =  HL;
        end

        if circularIndex(n+1,Nunit) == circularIndex(itN-1,Nunit)
        elseif circularIndex(n+1,Nunit) < circularIndex(itN-1,Nunit)
            for i = circularIndex(n+1,Nunit)+1:circularIndex(itN-1,Nunit)
                HL = updateLeft(HL,2,Al(i),[],[],Al(i));
            end
        else
            for i = circularIndex(n+1,Nunit)+1:circularIndex(itN-1,Nunit)+Nunit
                HL = updateLeft(HL,2,Al(circularIndex(i,Nunit)),[],[],Al(circularIndex(i,Nunit)));
            end
        end
        HLoutArr(n) = HL;
    end

    for n = 1:Nunit
        if circularIndex(n+1,Nunit) ~= itN
            HLinArr(n) = updateLeft(HLoutArr(n),2,Al(itN),[],[],Al(itN));
        end
    end

    HLout = sum(HLoutArr);
    HLin = sum(HLinArr);

    HRoutArr = QSpace(Nunit,1);
    HRinArr = HRoutArr;
    % Calculate 
    for n = 1:Nunit
        [HR, resNormHR(n,1)] =  gmres_HL_QSpace(Ar,C,h,n,nKrylov,tol,'right',true);
        % resNormHR(n,1)
        for iter = 1:iterNum-1
            [HR,resNormHR(n,iter+1)] = gmres_HL_QSpace(Ar,C,h,n,nKrylov,tol,'x0',HR,'right',true);
            % resNormHR(n,iter+1)
        end
        if n == itN
            HRinArr(n) =  HR;
        end

        if n == circularIndex(itN+1,Nunit)
        elseif n > circularIndex(itN+1,Nunit)
            for i = n-1:-1:circularIndex(itN+1,Nunit)
                HR = updateLeft(HR,2,permute(Ar(i),[2 1 3]),[],[],permute(Ar(i),[2 1 3]));
            end
        else
            for i = n-1 + Nunit:-1:circularIndex(itN+1,Nunit)
                HR = updateLeft(HR,2,permute(Ar(circularIndex(i,Nunit)),[2 1 3]),[],[],permute(Ar(circularIndex(i,Nunit)),[2 1 3]));
            end
        end
        HRoutArr(n) = HR;
    end

    for n = 1:Nunit
        if n ~= itN
            HRinArr(n) = updateLeft(HRoutArr(n),2,permute(Ar(itN),[2 1 3]),[],[],permute(Ar(itN),[2 1 3]));
        end
    end

    HRin = sum(HRinArr);
    HRout = sum(HRoutArr);
end
