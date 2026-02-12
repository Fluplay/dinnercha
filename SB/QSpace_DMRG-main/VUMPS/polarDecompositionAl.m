function [Al, errL] = polarDecompositionAl(AC, Cr)
    AcCd = contract(AC,2,Hconj(Cr),1);
    [U,~,V] = svdSB(AcCd,[1 2],[],1e-40);
    I = getIdentity(U,3,'-z');
    Al = contract(U,3,I,1);
    Al = contract(Al,3,V,1,[1 3 2]);
    errL = norm(AC-contract(Al,2,Cr,1,[1 3 2]));
end