function [M, S, Eeig] = centerDiagonalize(M,S,MPO,Hlr,Nkeep,varargin)
    % default vales of optional input parameters
    nKrylov = 5;
    % nKrylov = max(round(numel(M)/3),5);
    tol = 1e-8;
    
    % parse options
    while ~isempty(varargin)
        if numel(varargin) < 2
            error('ERR: Option should be set by a pair of option name and value.');
        end
        switch varargin{1}
            case 'Krylov'
                nKrylov = varargin{2};
                varargin(1:2) = [];
            case 'tol'
                tol = varargin{2};
                varargin(1:2) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end

    Nsite = numel(M);
    Ncenter = Nsite/2;
    Mtemp = contract(M(Ncenter),2,S,1,[1 3 2]);
    Aold = contract(Mtemp,2,M(Ncenter+1),1,[1,3,2,4]);
    [Anew,Eeig] = eigs_2site_GS_QSpace(Hlr(Ncenter),MPO(Ncenter),MPO(Ncenter+1),Hlr(Ncenter+3),Aold,nKrylov,tol);
    [M(Ncenter),S,M(Ncenter+1)] = svdMK(Anew,[1,3],Nkeep);
    M(Ncenter) = permute(M(Ncenter), [1 3 2]);
end