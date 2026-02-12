function EE = EEMPS(M)
%     <Describe>

%     Calculate entanglement entroy of every bond, including final bond (nonzero when MPS is mixed state)

%     <Input>
%     M : [1 x N QSpace array] MPS 

%     <Output>

%     EE : [1 x N] array of entanglement entropy, ith element is entanglement of bond M(i)--M(i+1)

%     Written by Subin Kim (01.02.2025)

    % M is left canaonical
    % output: 
    L = length(M);
    EE = zeros(1,L);
    [M,S] = canonFormMK(M,L,[]);
    S = S/norm(S);
    S2 = contract(S,2,S,'2*');
    EE(L) = SEntropySB(S2);
    % absorb S into M
    M(end) = legflip(M(end),2);
    M(end) = contract(M(end),2,S,1,[1 3 2]);
    for itN = (L:-1:2)
        itag = getitags(M(itN), 1);
        if ~isempty(itag) && itag(end) == '*'
            itag = itag(1:end-1);
        end
        [U, S, M(itN)] = svdMK(M(itN), 1, [], itag);
        M(itN-1) = contract(M(itN-1), 2, U*S, 1, [1 3 2]);
        S2 = contract(S,2,S,'2*');
        EE(itN-1) = SEntropy(S2);
        % EE(itN) = sum(Spart(~isnan(Spart)));
    end


end
