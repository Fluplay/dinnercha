function [Ar, errR] = polarDecompositionAr(AC, Cl)
    CdAC = contract(Hconj(Cl),2,AC,1);
    [U,~,V] = svdSB(CdAC,1,[],1e-40);
    I = getIdentity(U,2,'-z');
    Ar = contract(U,2,I,1);
    Ar = contract(Ar,2,V,1);
    errR = norm(AC-contract(Cl,2,Ar,1));
end