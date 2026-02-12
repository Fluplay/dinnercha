function [ACnew,Cl,Cr, resNormAC, resNormCl,resNormCr] = eigsSolver_QSpace(Al,Ar,C,h, HLout,HLin,HRin,HRout,itN,nKrylovEig,tol,iterEig)
    Nunit = numel(Al);

    % 초기화
    resNormAC = zeros(1, iterEig);
    resNormCl = zeros(1, iterEig);
    resNormCr = zeros(1, iterEig);
    
    ACold = contract(Al(itN),2,C(itN),1,[1,3,2]);
    [ACnew,~,resNormAC(1)] = eigs_HAC_QSpace(Al(circularIndex(itN-1,Nunit)),h(circularIndex(itN-1,Nunit)),Ar(circularIndex(itN+1,Nunit)),h(itN),HLout,HRout,ACold,nKrylovEig,tol);
    [Cl,~,resNormCl(1)] = eigs_HC_QSpace(Al(circularIndex(itN-1,Nunit)),h(circularIndex(itN-1,Nunit)),Ar(itN),HLout,HRin,C(circularIndex(itN-1,Nunit)),nKrylovEig,tol);
    [Cr,~,resNormCr(1)] = eigs_HC_QSpace(Al(itN),h(itN),Ar(circularIndex(itN+1,Nunit)),HLin,HRout,C(itN),nKrylovEig,tol);

    for iter = 1:iterEig-1
        [ACnew,~,resNormAC(iter+1)] = eigs_HAC_QSpace(Al(circularIndex(itN-1,Nunit)),h(circularIndex(itN-1,Nunit)),Ar(circularIndex(itN+1,Nunit)),h(itN),HLout,HRout,ACnew,nKrylovEig,tol);
        [Cl,~,resNormCl(iter+1)] = eigs_HC_QSpace(Al(circularIndex(itN-1,Nunit)),h(circularIndex(itN-1,Nunit)),Ar(itN),HLout,HRin,Cl,nKrylovEig,tol);
        [Cr,~,resNormCr(iter+1)] = eigs_HC_QSpace(Al(itN),h(itN),Ar(circularIndex(itN+1,Nunit)),HLin,HRout,Cr,nKrylovEig,tol);
    end
end