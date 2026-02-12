function Hlr = initialHlrCenter(M,S,Hs)
    % Hlr : environment Hamiltonian required for 2-site DMRG
    N = numel(M);
    Ncenter = N/2; % N is always even
    if N < 2
        error('ERR: chain is too short.');
    elseif N ~= numel(Hs)
        error('ERR: The lengths of ''M'' and ''Hs'' should be equal.');
    end
    
    if ~islegin(M(1), 1) || ~islegin(M(end), 2)
        error('please check the leg direction of the first or last MPS');
    end
    
    % Contractions of MPO and MPS tensors that represent the effective
    % Hamiltonian for the left/right parts of the chain: When the orthogonality
    % center (= central tensor in a site-canonical form) is at site n, then
    % Hlr(m+1) for m < n (m > n) is the left (right) part of the effective
    % Hamiltonian made of M(1:m) and Hs(1:m) [M(m:end) and Hs(m:end)].
    Hlr = QSpace(N+2,1);
    if isempty(M(1).info.itags{1})
        Hstart = getIdentity(M(1),1,Hs(1),3);
    else
        Hstart = getIdentity(M(1),1,Hs(1),3,strcat(M(1).info.itags{1},'*'));
    end
    Hstart = permute(Hstart,[3 1 2], 'conj'); % leg order: bottom (in), top (out), right (out)
    if isempty(M(end).info.itags{2})
        Hend = getIdentity(M(end),2,Hs(end),4, [1 3 2]); % leg order: bottom (in), top (out), left (in)
    else
        Hend = getIdentity(M(end),2,Hs(end),4,strcat(M(end).info.itags{2},'*'), [1 3 2]); % leg order: bottom (in), top (out), left (in)
    end
    Hlr(1) = Hstart; % leg order: bottom (in), top (out), right (out)
    Hlr(end) = Hend; % **NB!** leg order: bottom (in), top (out), left (in)
    
    for itN = (1:Ncenter-1)
        Hlr(itN+1) = updateLeft(Hlr(itN),3,M(itN),Hs(itN),4,M(itN));
    end
    for itN = (N-1:-1:Ncenter+1)
        Hlr(itN+2) = updateLeft(Hlr(itN+3), 3, permute(M(itN+1),[2 1 3]), permute(Hs(itN+1),[1 2 4 3]),4,permute(M(itN+1),[2 1 3]));
    end
end