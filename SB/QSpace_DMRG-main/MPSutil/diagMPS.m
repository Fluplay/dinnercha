function M = diagMPS(M)
    % <Decription>

    % Diagonalization of MPS M, if final bond dimension of M is not 1.

    % Written by Subin Kim (01.02.2025)

    Nsite = numel(M);
    N = normMPS(M,Nsite);
    [U,S,Ud] = svdQS(N,1,'stol',0);
    M(Nsite) = contract(M(Nsite),2,QSpace(Ud)',1,[1 3 2]);
    M(Nsite)
end