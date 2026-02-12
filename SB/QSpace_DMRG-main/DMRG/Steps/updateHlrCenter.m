function Hlr = updateHlrCenter(M,S,Hs,Hlr)
    % When iDMRG,so environment is not updated except about new MPS tensors.
    N = numel(M);
    Nsitebef = N-2;
    Ncenterbef = Nsitebef/2;
    Hlrbef = Hlr;
    if N < 2
        error('ERR: chain is too short.');
    elseif N ~= numel(Hs)
        error('ERR: The lengths of ''M'' and ''Hs'' should be equal.');
    end
    
    if ~islegin(M(1), 1) || ~islegin(M(end), 2)
        error('please check the leg direction of the first or last MPS');
    end
    
    Hlr = QSpace(N+2,1);


    Hlr(1:Ncenterbef) = Hlrbef(1:Ncenterbef);
    Hlr((Ncenterbef+5):(Nsitebef+4)) = Hlrbef((Ncenterbef+3):(Nsitebef+2));
    Hlr(Ncenterbef+1) = updateLeft(Hlr(Ncenterbef),3,M(Ncenterbef),Hs(Ncenterbef),4,M(Ncenterbef));
    Hlr(Ncenterbef+4) = updateLeft(Hlr(Ncenterbef+5), 3, permute(M(Ncenterbef+3),[2 1 3]), permute(Hs(Ncenterbef+3),[1 2 4 3]),4,permute(M(Ncenterbef+3),[2 1 3]));
end