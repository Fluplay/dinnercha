function [M,E0,Eiter,Sv,dw] = DMRG_GS_CBE_QSpace (M,Hs,Nkeep,Nsweep,delta,alpha,varargin)
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
    
    disptime('CBE-DMRG: ground state search');
    disptime(['# of sites = ',sprintf('%i',numel(Hs)), ...
        ', Nkeep = ',sprintf('%i',Nkeep),', # of sweeps = ',sprintf('%i',Nsweep),' x 2']);
    Sv = QSpace(Nsweep*2,N+1); % collection of singular value vectors
    [M,Sv(1,N)] = canonFormMK(M,N-1,[]); %Sv(1,n) is a bond between n-1 and n
    Sv(1,N) = diag(Sv(1,N));
    
    Eiter = zeros(N,2*Nsweep);
    dw = zeros(N-1,2*Nsweep);
        
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
    
    for itN = (1:N-1)
        Hlr(itN+1) = updateLeft(Hlr(itN),3,M(itN),Hs(itN),4,M(itN));
    end
    for itS = (1:Nsweep)
        % right -> left
        if itS > 1
            Sv(2*itS-1,N) = Sv(2*itS-2,N);
            if alpha > 1
                Nkeep = ceil(Nkeep*alpha);
                fprintf('Nkeep = %d\n',Nkeep);
            end
        end

        for itN = (N:-1:2)
            Atr = CBE_shrewd(Hlr(itN-1),Hlr(itN+2),M(itN-1),M(itN),Hs(itN-1),Hs(itN),Sv(2*itS-1,itN),Nkeep,delta);
            Aex=cat(M(itN-1),Atr,2);
            Hex = updateLeft(Hlr(itN-1),3,Aex,Hs(itN-1),4,Aex);
            Ainit = contract(conj(Aex),[1 3],M(itN-1),[1 3]);
            Ainit = contract({Ainit,2,diag(Sv(2*itS-1,itN)),1},2,M(itN),1);
            [M(itN),Eiter(N+1-itN,2*itS-1)] = eigs_1site_GS_QSpace(Hex,Hs(itN),Hlr(itN+2),Ainit,nKrylov,tol);
            % SVD
            [U,Sv(2*itS-1, itN),M(itN),dw(itN-1,2*itS-1)] = svdSB(M(itN),1,Nkeep,0,M(itN).info.itags{1});
            Sv(2*itS-1, itN) = diag(Sv(2*itS-1, itN));
            %normalization
            % Sv(2*itS-1, itN) = Sv(2*itS-1, itN)/norm(Sv(2*itS-1, itN));

            % update the next tensor
            if itN > 2
                M(itN-1)=contract(Aex,2,U*diag(Sv(2*itS-1,itN)),1,[1 3 2]);
                [U,Sv(2*itS-1,itN-1),M(itN-1)] =svdSB(M(itN-1),1,Nkeep,0,M(itN-1).info.itags{1});
                Sv(2*itS-1,itN-1) = diag(Sv(2*itS-1,itN-1));
                M(itN-2) = contract(M(itN-2),2,U,1,[1 3 2]);
            else
                M(itN-1) = contract(Aex,2,U,1,[1 3 2]);
            end
            Hlr(itN+1) = updateLeft(Hlr(itN+2), 3, permute(M(itN),[2 1 3]), permute(Hs(itN),[1 2 4 3]),4,permute(M(itN),[2 1 3]));
        end
        disptime(['Sweep #',sprintf('%i/%i',[2*itS-1,2*Nsweep]),' (right -> left): Energy = ', ...
            sprintf('%.16g',Eiter(N-1,2*itS-1))]);

        Sv(2*itS,2) = Sv(2*itS-1,2);
        % left -> right
        for itN = (1:N-1)
            % Atr = CBE_shrewd(Hlr(itN+3),Hlr(itN),permute(M(itN+1),[2 1 3]),permute(M(itN),[2 1 3]),permute(Hs(itN+1),[1 2 4 3]),permute(Hs(itN),[1 2 4 s]),Sv(2*itS,itN+1),Nkeep,delta);
            Atr = CBE_shrewd(Hlr(itN+3),Hlr(itN),permute(M(itN+1),[2 1 3]),permute(M(itN),[2 1 3]),permute(Hs(itN+1),[1 2 4 3]),permute(Hs(itN),[1 2 4 3]),permute(Sv(2*itS,itN+1),[2 1]),Nkeep,delta);

            Atr = permute(Atr,[2 1 3]);
            Aex = cat(M(itN+1),Atr,1);
            Hex = updateLeft(Hlr(itN+3),3,permute(Aex,[2 1 3]),permute(Hs(itN+1),[1 2 4 3]),4,permute(Aex,[2 1 3]));
            Ainit = contract(M(itN+1),[2 3],conj(Aex),[2 3]);
            Ainit = contract(M(itN),2,diag(Sv(2*itS,itN+1))*Ainit,1,[1 3 2]);
            [M(itN), Eiter(itN,2*itS)] = eigs_1site_GS_QSpace(Hlr(itN),Hs(itN),Hex,Ainit,nKrylov,tol);
            [M(itN),Sv(2*itS,itN+1),Vd] = svdSB(M(itN),[1 3],Nkeep,0,M(itN).info.itags{2});
            Sv(2*itS,itN+1) = diag(Sv(2*itS,itN+1));
            M(itN) = permute(M(itN),[1 3 2]);
            if itN < N-1
                M(itN+1) = contract({diag(Sv(2*itS,itN+1)),'!1',Vd},'!1',Aex);
                [M(itN+1),Sv(2*itS,itN+2),Vd] = svdSB(M(itN+1),[1 3],Nkeep,0,M(itN+1).info.itags{2});
                Sv(2*itS,itN+2) = diag(Sv(2*itS,itN+2));
                M(itN+1) = permute(M(itN+1),[1 3 2]);
                M(itN+2) = contract(Vd,2,M(itN+2),1);
            else
                % Sv(N) should be 1, as the MPS norm; To keep the arrow direction of the first leg of the leftmost MPS to in, absorb Sv and Vd into M(1)
                M(itN+1) = contract(Vd,2,Aex,1);
                % M(itN+1) = contract(M(itN+1),2,Sv(2*itS,itN+1),1);
            end
            Hlr(itN+1) = updateLeft(Hlr(itN), 3, M(itN), Hs(itN), 4, M(itN));
        end
        disptime(['Sweep #',sprintf('%i/%i',[2*itS,2*Nsweep]),' (left -> right) : Energy = ', ...
            sprintf('%.16g',Eiter(N-1,2*itS))]);
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

    function Atr = CBE_shrewd(HL,HR,Al,Ar,Hl,Hr,Sb,Nkeep,delta)
        % bond dimension numbers
        D = dim(Sb,2);
        w = dim(Hl,4);
        Dp = ceil(D/w);
        Dt = ceil(Nkeep*delta);
        % Al;
        % contract(conj(Al),2,Al,2)
        % id1 = getIdentity(Al,1);
        % id2 = getIdentity(Al,2);
        % id = directprod(id1,id2);
        % id = permute(id,[1 3 2 4])
        % contract HL, Al and Hl
        TL = contract(HL,2,Al,1);
        TL = contract(TL,[2 4],Hl,[3 2],[1 3 2 4]);

        % contract TL with the projector made of Al
        TLk = contract(conj(Al),[1 3],TL,[1 2]);
        TLk = contract(Al,2,TLk,1);

        % leg order: left of Al, bottom of Al, right of AL, right of Hl
        
        % subtraction between TL and TLk, as the left half of the first diagram at the upper-left corner
        TLd = TL - TLk;
        % contract HR, Ar, and Hr
        TR = contract(HR,2,Ar,2);
        TR = contract(TR,[2 4],Hr,[4 2],[1 3 2 4]);
        % leg order : bottom of HR, bottom of Hr, left of Ar, left of Hl

        % contract HR with the projector made of Ar
        TRk = contract(conj(Ar),[2 3],TR,[1 2]);
        TRk = contract(Ar,1,TRk,1);

        % substraction
        TRd = TR - TRk;
        % SVD for the part shaded in pink
        T = contract(diag(Sb),2,TRd,3);
        [U,S] = svdMK(T,1,D);
        % SVD for the part shaded in blue
        T = contract(TLd, 3, contract(U,2,S,1),1);
        [U,S] = svdSB(T,[1 2 3],Dp,1e-10);
        if isempty(S.data)
            Atr =  QSpace();
            return
        end
        % SVD for the part shaded in red
        T = contract(U,4,S,1);
        [Apr,S] = svdSB(T,[1 2],Dp*w,1e-10); %Dp*w
        if isempty(S.data)
            Atr =  QSpace();
            return
        end
        % SVD for the part shaded in oragne
        T = contract(conj(Apr),[1 2],TL,[1 2]);
        T = contract(T,2,diag(Sb),1);
        T = contract(T,[2 3],TRd,[4 3]);
        if isempty(T.data)
            Atr =  QSpace();
            return
        end
        U = svdSB(T,1,Dt,1e-10);

        Atr = contract(Apr,3,U,1,[1 3 2]);
        
    end