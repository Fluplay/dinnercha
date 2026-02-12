function [M,E0,Eiter,Sv,Hlr] = DMRG_GS_1site_QSpace (M,Hs,Nkeep,Nsweep,varargin)
% < Description >
%
% [M,E0,Eiter,Sv] = DMRG_GS_1site (M,Hs,Nkeep,Nsweep [,'Krylov',nKrylov] [,'tol',tol])
%
% Single-site density-matrix renormalization group (DMRG) calculation to
% search for the ground state and its energy of a one-dimensional system,
% whose Hamiltonian is given by the matrix product operator.
%
% < Input >
% M : [1 x N QSpace array] MPS as the initial guess of
%       the variational search. Each M(n) is a rank-3 tensor acting on site
%       n. Their legs are ordered as left-right-bottom(=physical). The
%       length M, i.e., numel(Hs), defines the system length N.
%       **NB!**
%       The direction of the first leg of the MPS(1) and the second leg of
%       the MPS(end) must be in.
% Hs : [1 x N QSpace array] MPO of the Hamiltonian.
%       Each Hs(n) is a rank-4 tensor acting on site n. The order of legs
%       of Hs(n) is bottom-top-left-right, where the bottom (top) leg
%       contracts to the physical leg of bra (ket) tensor. 
% Nkeep : [numeric] Maximum bond dimension of the MPS to keep.
% Nsweep : [numeric] Number of round trips of sweeping. There will be
%       Nsweep pairs of right-to-left sweep and left-to-right sweep. That
%       is, the first sweep is right-to-left and the last is left-to-right.
%
% < Option >
% 'Krylov', .. : [numeric] The maximum dimension of the Krylov subspace to
%       be considered, within the Lanczos method used for updating each MPS
%       tensor. That is, a tridiagonal matrix of size n-by-n (at maximal)
%       is diagonalized in each tensor update, with n set by 'Krylov',.. .
%       (Default: 5)
% 'tol', .. : [numeric] The tolerance for elements on the +-1 diagonals
%       (those next to the main diagonal) within the Lanczos method. If an
%       element is smaller than the tolerance, then the Lanczos
%       tridiagonalization procedure stops and only the part of the
%       tridiagonal matrix constructed so far is diagonalized.
%       (Default: 1e-8)
%
% < Output >
% M : [1 x N QSpace array] The result MPS which is obtained variationally to
%       have the minimum expectation value of the Hamiltonian Hs. It is in
%       a left-canonical form, since the last sweep is left-to-right.
% E0 : [numeric] The energy of M.
% Eiter : [N x (2*Nsweep) numeric array] Each element Eiter(m,n) means the
%       variational energy in the m-th iteration within the n-th sweep.
%       Odd n is for right-to-left sweep and even n for left-to-right
%       sweep. Note that the iteration index m matches with the site index
%       for left-to-right sweep; the iteration m corresponds to the site
%       (N+1-m) for right-to-left sweep.
% Sv : [1 x (N+1) QSpace array] Sv(n) contains a vector of singular values on
%       the bond between sites n-1 and n, for 1 < n < N. Sv(1) and Sv(end)
%       are the norms of the MPS, sitting on the left- and rightmost legs
%       (that are dummy), respectively.
%
% Written by Minsoo Kim (Dec. 06,2022)


tobj = tic2;

% default vales of optional input parameters
nKrylov = 5;
tol = 1e-8;
canonical = false;
makeHlr = true;

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
        case 'Canonical'
            canonical = varargin{2};
            varargin(1:2) = [];
        case 'Hlr'
            Hlr = varargin{2};
            makeHlr = false;
            varargin(1:2) = [];  
        otherwise
            error('ERR: Unknown input.');
    end
end

% % sanity check for input and option
N = numel(M);
if N < 2
    error('ERR: chain is too short.');
elseif N ~= numel(Hs)
    error('ERR: The lengths of ''M'' and ''Hs'' should be equal.');
end

if ~islegin(M(1), 1) || ~islegin(M(end), 2)
    error('please check the leg direction of the first or last MPS');
end

disptime('Single-site DMRG: ground state search');
disptime(['# of sites = ',sprintf('%i',numel(Hs)), ...
    ', Nkeep = ',sprintf('%i',Nkeep),', # of sweeps = ',sprintf('%i',Nsweep),' x 2']);

if ~canonical
    M = canonFormMK(M,N,[]);
end

Eiter = zeros(N,2*Nsweep);

Sv = QSpace(Nsweep*2,N+1); % collection of singular value vectors

% Contractions of MPO and MPS tensors that represent the effective
% Hamiltonian for the left/right parts of the chain: When the orthogonality
% center (= central tensor in a site-canonical form) is at site n, then
% Hlr(m+1) for m < n (m > n) is the left (right) part of the effective
% Hamiltonian made of M(1:m) and Hs(1:m) [M(m:end) and Hs(m:end)].
if makeHlr
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

    for itN = (1:N)
        Hlr(itN+1) = updateLeft(Hlr(itN),3,M(itN),Hs(itN),4,M(itN));
    end

end

for itS = (1:Nsweep)
    % right -> left
    for itN = (N:-1:1)
        [M(itN),Eiter(N+1-itN,2*itS-1)] = eigs_1site_GS_QSpace(Hlr(itN),Hs(itN),Hlr(itN+2),M(itN),nKrylov,tol);
        [M(itN),Sv(2*itS-1, itN),Vd] = svdQS(M(itN),1,'Nkeep',-1,'stol',0);
        M(itN)=QSpace(M(itN)); Sv(2*itS-1, itN)=QSpace(Sv(2*itS-1, itN)); Vd=QSpace(Vd); 

        % update the next tensor
        if itN > 1
            M(itN-1)=contract(M(itN-1),2,Vd,2);
            M(itN-1)=contract(M(itN-1),3,diag(Sv(2*itS-1,itN)),2, [1 3 2]);
        else
            % Sv(N+1) should be 1, as the MPS norm; To keep the arrow direction of the second leg of the rightmost MPS to in, absorb Sv and Vd into M(N)
            M(itN) = contract(Vd, 1, contract(diag(Sv(2*itS-1, itN)), 1, M(itN), 1), 1);
        end
        Hlr(itN+1) = updateLeft(Hlr(itN+2), 3, permute(M(itN),[2 1 3]), permute(Hs(itN),[1 2 4 3]),4,permute(M(itN),[2 1 3]));
    end
    disptime(['Sweep #',sprintf('%i/%i',[2*itS-1,2*Nsweep]),' (right -> left) : Energy = ', ...
        sprintf('%.16g',Eiter(N,2*itS-1))]);
    
    % left -> right
    for itN = (1:N)
        [M(itN), Eiter(itN,2*itS)] = eigs_1site_GS_QSpace(Hlr(itN),Hs(itN),Hlr(itN+2),M(itN),nKrylov,tol);
        [M(itN),Sv(2*itS,itN+1),Vd] = svdQS(M(itN),2,'Nkeep',-1,'stol',0);
        M(itN) = QSpace(M(itN)); Sv(2*itS,itN+1) = QSpace(Sv(2*itS, itN+1)); Vd = QSpace(Vd);
        
        if itN < N
            M(itN+1) = contract({diag(Sv(2*itS,itN+1)),'!1',Vd},'!1',M(itN+1));
        else
            % Sv(1) should be 1, as the MPS norm; To keep the arrow direction of the first leg of the leftmost MPS to in, absorb Sv and Vd into M(1)
            M(itN) =contract(contract(M(itN), 2, diag(Sv(2*itS,itN+1)), 1), 3, Vd, 1, [1 3 2]);
        end
        Hlr(itN+1) = updateLeft(Hlr(itN), 3, M(itN), Hs(itN), 4, M(itN));
    end
    disptime(['Sweep #',sprintf('%i/%i',[2*itS,2*Nsweep]),' (left -> right) : Energy = ', ...
        sprintf('%.16g',Eiter(N,2*itS))]);
end

E0 = Eiter(N,2*itS); % take the last value
    
toc2(tobj,'-v');

end

function [Anew, Enew] = eigs_1site_GS_QSpace (Hleft,Hcen,Hright,Aold,nKrylov,tol)
% < Description >
%
% [Anew, Enew] = eigs_1site_GS_QSpace (Hleft,Hcen,Hright,Aold,nKrylov,tol)
%
% Update an MPS tensor acting on a single site, by solving the effective
% Hamiltonian via the Lanczos method.
%
% < Input >
% Hleft : [rank-3 QSpace tensor] Left part of the effective Hamiltonian. Its legs
%       are ordered as bottom-top-right.
% Hcen : [rank-4 QSpace tensor] Center part of the effective Hamiltonian. It is
%       indeed an MPO tensor for the current site.
% Hright : [rank-3 QSpace tensor] Right part of the effective Hamiltonian. Its
%       legs are ordered as bottom-top-left.
% Aold : [rank-3 QSpace tensor] Current ket tensor.
% nKrylov, tol : Options for the Lacnzos method. See the documentation for
%       the parent function for details.
%
% < Output >
% Anew : [rank-3 QSpace tensor] The approximate ground state of the effective
%       Hamiltonian obtained by the Lanczos method, using the Krylov
%       subspace with dimension up to "nKrylov".
% Enew : [numeric] Expectation value of the effective Hamiltonian with
%       respect to "Anew".

As = QSpace(nKrylov,1);

As(1) = Aold/norm(Aold);

% normalize; insert "abs" to avoid numerical noise

alphas = zeros(nKrylov,1); % main diagonal elements
betas  = zeros(nKrylov-1,1); % +-1 diagonal elements
cnt = 0; % counter for Krylov subspace dimension

for itn = (1:nKrylov)
    % "matrix-vector" multiplication. Contraction order can affect to the caculation time!
    Amul = contract(Hleft, 2, As(itn), 1);
    Amul = contract(Amul, [2 4], Hcen, [3 2]);
    Amul = contract(Amul, [2 4], Hright, [2 3], [1 3 2]);
    alphas(itn) = real(getscalarMK(contract(conj(As(itn)),(1:3),Amul,(1:3))));

    cnt = cnt+1;
    if itn < nKrylov
        for it2 = (1:2) % do twice to further reduce numerical noise
            for itK = 1:itn
                T = contract(conj(As(itK)), [1 2 3] ,Amul, [1 2 3]);
                T = getscalarMK(T)*As(itK);
                Amul = Amul - T;
            end
        end

        Anorm = norm(Amul);
        if Anorm < tol % for numerical stability
            break;
        end
        As(itn+1) = Amul/Anorm; % normalize
        betas(itn) = Anorm;
    end
    % HKrylovx = diag(betas(1:cnt-1),-1);
    % HKrylovx = HKrylovx + HKrylovx' + diag(alphas(1:cnt));
    % [V,D] = eig(HKrylovx);
    % D = diag(D);
    % D(1:min(cnt,5))
    % betas(min(cnt,nKrylov-1))

end

Hkrylov = diag(betas(1:cnt-1),-1);
Hkrylov = Hkrylov + Hkrylov' + diag(alphas(1:cnt));

[V,D] = eig(Hkrylov);
[~,minid] = min(diag(D));
Anew=V(1,minid)*As(1);
for itK = 2:cnt
    Anew=Anew+V(itK,minid)*As(itK);
end
% AA = contract(conj(Anew),[1 3],Anew,[1 3])
% for i = 1:length(AA.data)
%     AA.data{i}
% end
% compute the epectation value of the effective Hamiltonian with respect to "Anew"

Enew = contract(Hleft, 2, Anew, 1);
Enew = contract(Enew, [2 4], Hcen, [3 2]);
Enew = contract(Enew, [2 4], Hright, [2 3], [1 3 2]);
Enew = contract(conj(Anew), [1 2 3], Enew, [1 2 3]);
Enew = real(getscalarMK(Enew));
end