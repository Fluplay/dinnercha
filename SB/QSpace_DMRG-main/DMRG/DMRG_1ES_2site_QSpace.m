function [M,E0,Eiter,Sv] = DMRG_1ES_2site_QSpace (M,M0, Hs,Nkeep,Nsweep,varargin)
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
    % Updated by Subin Kim(Jul. 15. 2024)
    
    
    tobj = tic2;
    
    % default vales of optional input parameters
    nKrylov = 5;
    tol = 1e-8;
    near0=1e-40;
    
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
    
    disptime('Double-site DMRG: 1st excited state search');
    disptime(['# of sites = ',sprintf('%i',numel(Hs)), ...
        ', Nkeep = ',sprintf('%i',Nkeep),', # of sweeps = ',sprintf('%i',Nsweep),' x 2']);
    
    % M = canonFormMK(M,N,[]);
    % M0 = canonFormMK(M0,N,[]);
    Eiter = zeros(N-1,2*Nsweep);
    
    Sv = QSpace(Nsweep*2,N+1); % collection of singular value vectors
    
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


    % Overlap between MPS and ground M0.
    Olr = QSpace(N+2,1);
    if isempty(M(1).info.itags{1})
        % Ostart = getIdentity(M(1),1,Hs(1),3);
        Ostart = getIdentity(M0(1),1);
    else
        % Ostart = getIdentity(M(1),1,Hs(1),3,strcat(M(1).info.itags{1},'*'));
        % Ostart = getIdentity(M0(1),1,strcat(M0(1).info.itags{1},'*'));
        A = getIdentity(M(1),1,legflip(M0(1),1),1,strcat(M0(1).info.itags{1}));
        A = legflip(A,2);
        v = getvac(A,'-1d');
        A = getsub(A,v.Q{1},3);
        for itM =(1:numel(A.data))
            A.data{itM} = A.data{itM}(:,:,1);
        end
        A = contract(A,3,v,1);
        Ostart = A;

    end
    % Ostart = permute(Ostart,[3 1 2], 'conj'); % leg order: bottom (in), top (out), right (out)
    
    if isempty(M(end).info.itags{2})
        % Oend = getIdentity(M0(end),2); % leg order: bottom (in), top (out), left (in)
    else
        Oend = getIdentity(M0(end),2,strcat(M0(end).info.itags{2})); % leg order: bottom (in), top (out), left (in)
        % M(end)
        % M0(end)
        % A = getIdentity(M(end),2,legflip(M0(end),2),2,strcat(M0(end).info.itags{2}));
        % A = legflip(A,2);
        % v = getvac(A,'-1d');
        % A = getsub(A,v.Q{1},3);
        % for itM =(1:numel(A.data))
        %     % % itM
        %     % D = A.data{itM};
        %     % V = [];
        %     % for i = 1:size(D,3)
        %     %     if i == 1
        %     %         V = D(:,:,i);
        %     %     else
        %     %         V = V+ D(:,:,i);
        %     %     end
        %     % end
        %     A.data{itM} = A.data{itM}(:,:,1);
        %     A.data{itM}
        %     % V
        % end
        % A = contract(A,3,v,1);
        % Oend = A;
        % Oend
        % Oend.data{1}
        % Oend.data{2}
    end
    % Ostart
    Olr(1) = Ostart; % leg order: bottom (in), top (out), right (out)
    Olr(end) = Oend; % **NB!** leg order: bottom (in), top (out), left (in)
    for itN = (1:N)
        if itN == 1
            % Hlr(itN+1) = updateLeft([],[],M(itN),permute(Hs(itN),[1 2 4 3]),3,M(itN));
            Hlr(itN+1) = updateLeft(Hlr(itN),3,M(itN),Hs(itN),4,M(itN));
        else
            Hlr(itN+1) = updateLeft(Hlr(itN),3,M(itN),Hs(itN),4,M(itN));
        end
        Olr(itN+1) = updateLeft(Olr(itN),2,M(itN),[],[],M0(itN));
        % if norm(Olr(itN+1)) == 0 % care of empty QSpace case
        %     Olr(itN+1) = conj(getIdentity(M(itN) + M0(itN),2)*near0);
        % end
        
    end
    for itS = (1:Nsweep)
        % right -> left
        for itN = (N-1:-1:1)
            % if itN>=7
                % M0(itN)
                % M0(itN+1)
                % Olr(itN)
                % Olr(itN+3)
                % M(itN)
                % M(itN+1)
            % end
            % itN
            % Olr(itN)
            Aorth = contract(M0(itN),2,M0(itN+1),1,[1 3 2 4]); 
            Aorth = contract(Olr(itN), 2, Aorth, 1) ;
            Aorth = contract(Aorth,2,Olr(itN+3),2,[1 4 2 3]);
            % take care of itag vanishing problem of contraction
            midtag = M(itN).info.itags{2};
    
            Aold = contract(M(itN),2,M(itN+1),1,[1,3,2,4]);
            % Aold
            [Anew,Eiter(N-itN,2*itS-1)] = eigs_2site_1ES_QSpace(Hlr(itN),Hs(itN),Hs(itN+1),Hlr(itN+3),Aold,Aorth,nKrylov,tol);
            % Anew
            % SVD; now we can safely truncate small singular values
            [M(itN),Sv(2*itS-1, itN+1),M(itN+1)] = svdMK(Anew,[1,3],Nkeep);
            M(itN) = permute(M(itN), [1 3 2]);
    
    
            % care tags
            M(itN) = setitags(M(itN),2,midtag);
            Sv(2*itS-1, itN+1) = setitags(Sv(2*itS-1, itN+1), [1 2], {midtag, midtag});
            M(itN+1) = setitags(M(itN+1),1,midtag);

    
            % update the next tensor
            M(itN) = contract(M(itN),2,Sv(2*itS-1,itN+1),1,[1,3,2]); 
    
            % update the Hamiltonian in effective basis
            Hlr(itN+2) = updateLeft(Hlr(itN+3), 3, permute(M(itN+1),[2 1 3]), permute(Hs(itN+1),[1 2 4 3]),4,permute(M(itN+1),[2 1 3]));
            % update overlap in effective basis
            % AorthMid = contract(M0(itN),2,M0(itN+1),1,[1 3 2 4]); 
            % M0(itN+1) = legflip(M0(itN+1),1);
            % [M0(itN),S0,M0(itN+1)] = svdMK(AorthMid,[1,3],[]);
            % M0(itN) = permute(M0(itN), [1 3 2]);
            % M0(itN) = setitags(M0(itN),2,midtag);
            % S0 = setitags(S0, [1 2], {midtag, midtag});
            % M0(itN+1) = setitags(M0(itN+1),1,midtag);
            % M0(itN) = contract(M0(itN),2,S0,1,[1,3,2]);
            M0(itN+1) = legflip(M0(itN+1),1);
            M0(itN) = legflip(M0(itN),2);
            Olr(itN+2) = updateLeft(Olr(itN+3),2,permute(M(itN+1),[2 1 3]),[],[],permute(M0(itN+1),[2 1 3]));
            % if norm(Olr(itN+2)) == 0 %care of empty QSpace case
            %     Olr(itN+2) = getIdentity(M(itN+1) + M0(itN+1),1)*near0;
            % end        
        end
        disptime(['Sweep #',sprintf('%i/%i',[2*itS-1,2*Nsweep]),' (right -> left) : Energy = ', ...
            sprintf('%.16g',Eiter(N-1,2*itS-1))]);
        % left -> right
        for itN = (1:N-1)
            % itN
            % Olr(itN+3)
            Aorth = contract(M0(itN),2,M0(itN+1),1,[1 3 2 4]); 
            Aorth = contract(Olr(itN), 2, Aorth, 1); 
            Aorth = contract(Aorth,2,Olr(itN+3),2,[1 4 2 3]);

            midtag = M(itN).info.itags{2};
            Aold = contract(M(itN),2,M(itN+1),1,[1,3,2,4]);
            % if itN == 7
            %     contract(conj(M(itN+1)),[1 3],M(itN+1),[1 3]).data{1}
            %     contract(conj(Aold), [1 3 4],Aold, [1 3 4]).data{1}
            % end
            [Anew, Eiter(itN,2*itS)] = eigs_2site_1ES_QSpace(Hlr(itN),Hs(itN),Hs(itN+1),Hlr(itN+3),Aold,Aorth, nKrylov,tol);
    
            % SVD; now we can safely truncate smal singular values
    
            [M(itN),Sv(2*itS,itN+1),M(itN+1)] = svdMK(Anew,[1 3],Nkeep);
            M(itN) = permute(M(itN),[1 3 2]);   
            % midtag inserting
            M(itN) = setitags(M(itN),2,midtag);
            Sv(2*itS, itN+1) = setitags(Sv(2*itS, itN+1), [1 2], {midtag, midtag});
            M(itN+1) = setitags(M(itN+1),1,midtag);
            % update the next tensor
            M(itN+1) = contract(Sv(2*itS,itN+1),2,M(itN+1),1);
            % update Hlr
            Hlr(itN+1) = updateLeft(Hlr(itN), 3, M(itN), Hs(itN), 4, M(itN));
            % update overlap
            % AorthMid = contract(M0(itN),2,M0(itN+1),1,[1 3 2 4]); 
            % [M0(itN),S0,M0(itN+1)] = svdMK(AorthMid,[1,3],[]);
            % M0(itN) = permute(M0(itN), [1 3 2]);
            % M0(itN) = setitags(M0(itN),2,midtag);
            % S0 = setitags(S0, [1 2], {midtag, midtag});
            % M0(itN+1) = setitags(M0(itN+1),1,midtag);
            % M0(itN+1) = contract(S0,2,M0(itN+1),1);
            M0(itN+1) = legflip(M0(itN+1),1);
            M0(itN) = legflip(M0(itN),2);
            Olr(itN+1) = updateLeft(Olr(itN),2,M(itN),[],[],M0(itN));
        end
        disptime(['Sweep #',sprintf('%i/%i',[2*itS,2*Nsweep]),' (left -> right) : Energy = ', ...
            sprintf('%.16g',Eiter(N-1,2*itS))]);
    end
    
    E0 = Eiter(N-1,2*itS); % take the last value
        
    toc2(tobj,'-v');
    
    end
    
    function [Anew, Enew] = eigs_2site_1ES_QSpace (Hleft,Hcen1,Hcen2,Hright,Aold,Aorth,nKrylov,tol)
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
    % contract(conj(Aorth),[1 2 3 4],Aold,[1 2 3 4])
    As = QSpace(nKrylov+1,1);
    Aorth = Aorth/sqrt(real(getscalarMK(contract(conj(Aorth),[1 2 3 4],Aorth,[1 2 3 4]))));
    Aold = Aold-Aorth*getscalarMK(contract(conj(Aorth),[1 2 3 4],Aold,[1 2 3 4]));
    Aold = Aold-Aorth*getscalarMK(contract(conj(Aorth),[1 2 3 4],Aold,[1 2 3 4]));
    As(1) = Aorth;
    As(2) = Aold/sqrt(real(getscalarMK(contract(conj(Aold),[1 2 3 4],Aold,[1 2 3 4])))) ;
    % As = QSpace(nKrylov+1,1);
    % Aorth = Aorth/norm(Aorth);
    % Aold = Aold-Aorth*norm(contract(conj(Aorth),[1 3 4],Aold,[1 3 4]));
    % Aold = Aold-Aorth*norm(contract(conj(Aorth),[1 3 4],Aold,[1 3 4]));
    % As(1) = Aorth;
    % As(2) = Aold/norm(Aold);
    % norm(contract(conj(Aorth),[1 3 4],Aold,[1 3 4]))
    
    alphas = zeros(nKrylov,1); % main diagonal elements
    betas  = zeros(nKrylov-1,1); % +-1 diagonal elements
    cnt = 0; % counter for Krylov subspace dimension
    
    for itn = (1:nKrylov)
        % "matrix-vector" multiplication. Contraction order can affect to the caculation time!
        Amul = contract(Hleft, 2, As(itn+1), 1);
        Amul = contract(Amul, [2 4], Hcen1, [3 2]);
        Amul = contract(Amul, [3 5], Hcen2, [2 3]);
        Amul = contract(Amul, [2 5], Hright, [2 3], [1 4 2 3]);
        alphas(itn) = real(getscalarMK(contract(conj(As(itn+1)),[1 2 3 4],Amul,[1 2 3 4])));
        % alphas(itn) = norm(contract(conj(As(itn+1)),[1 3 4],Amul,[1 3 4]));
    
        cnt = cnt+1;
        if itn < nKrylov
            % orthogonalize, to get the next Krylov vector
            for it2 = (1:2) % do twice to further reduce numerical noise
                for itK = 1:itn+1
                    T = contract(conj(As(itK)), [1 2 3 4] ,Amul, [1 2 3 4]);
                    T = getscalarMK(T)*As(itK);
                    % T = contract(conj(As(itK)), [1 3 4] ,Amul, [1 3 4]);
                    % T = norm(T)*As(itK);
                    Amul = Amul - T;
                end
            end
    
            Anorm = sqrt(real(getscalarMK(contract(conj(Amul),[1 2 3 4],Amul,[1 2 3 4]))));
            % Anorm = norm(Amul);
            if Anorm < tol % for numerical stability
                break;
            end
            As(itn+2) = Amul/Anorm; % normalize
            betas(itn) = Anorm;
        end
    end
    
    Hkrylov = diag(betas(1:cnt-1),-1);
    Hkrylov = Hkrylov + Hkrylov' + diag(alphas(1:cnt));
    
    [V,D] = eig(Hkrylov);
    [~,minid] = min(diag(D));
    Anew=V(1,minid)*As(2);
    for itK = 3:cnt+1
        Anew=Anew+V(itK-1,minid)*As(itK);
    end
    % Aorth
    % Anew
    % contract(conj(Aorth),[1 3 4],Aorth,[1 3 4]).data{1}
    % contract(conj(Anew),[1 3 4],Anew,[1 3 4]).data{1}
    % try
    % contract(conj(Aorth),[1 2 3 4],Anew,[1 2 3 4])
    % catch
    % end
    % norm(contract(conj(Aorth),[1 3 4],Anew,[1 3 4]))
    % compute the epectation value of the effective Hamiltonian with respect to "Anew"
    Amul = contract(Hleft, 2, Anew, 1);
    Amul = contract(Amul, [2 4], Hcen1, [3 2]);
    Amul = contract(Amul, [3 5], Hcen2, [2 3]);
    Amul = contract(Amul, [2 5], Hright, [2 3], [1 4 2 3]);
    Amul = contract(conj(Anew), [1 2 3 4], Amul, [1 2 3 4]);
    Enew = real(getscalarMK(Amul));

    
    
    
    end