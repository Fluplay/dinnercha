function [Anew, Enew] = eigs_2site_GS_QSpace (Hleft,Hcen1,Hcen2,Hright,Aold,nKrylov,tol)
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
    
    As(1) = Aold/sqrt(real(getscalarMK(contract(conj(Aold),[1 2 3 4],Aold,[1 2 3 4]))));
    
    % normalize; insert "abs" to avoid numerical noise
    
    alphas = zeros(nKrylov,1); % main diagonal elements
    betas  = zeros(nKrylov-1,1); % +-1 diagonal elements
    cnt = 0; % counter for Krylov subspace dimension
    
    for itn = (1:nKrylov)
        % "matrix-vector" multiplication. Contraction order can affect to the caculation time!
        Amul = contract(Hleft, 2, As(itn), 1);
        Amul = contract(Amul, [2 4], Hcen1, [3 2]);
        Amul = contract(Amul, [3 5], Hcen2, [2 3]);
        Amul = contract(Amul, [2 5], Hright, [2 3], [1 4 2 3]);
        alphas(itn) = real(getscalarMK(contract(conj(As(itn)),[1 2 3 4],Amul,[1 2 3 4])));
        % alphas(itn) = real(getscalarMK(contract(conj(As(itn)),(1:4),Amul,(1:4))));
    
        cnt = cnt+1;
        if itn < nKrylov
            % orthogonalize, to get the next Krylov vector
            for it2 = (1:2) % do twice to further reduce numerical noise
                for itK = 1:itn
                    T = contract(conj(As(itK)), [1 2 3 4] ,Amul, [1 2 3 4]);
                    T = getscalarMK(T)*As(itK);
                    Amul = Amul - T;
                end
            end
            Anorm = sqrt(real(getscalarMK(contract(conj(Amul),[1 2 3 4],Amul,[1 2 3 4]))));
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
    % diag(D)
    [~,minid] = min(diag(D));
    Anew=V(1,minid)*As(1);
    for itK = 2:cnt
        Anew=Anew+V(itK,minid)*As(itK);
    end
    % AA = contract(conj(Anew),[1 3 4],Anew,[1 3 4])
    % for i = 1:length(AA.data)
    %     AA.data{i}
    % end
    % compute the epectation value of the effective Hamiltonian with respect to "Anew"
    Amul = contract(Hleft, 2, Anew, 1);
    Amul = contract(Amul, [2 4], Hcen1, [3 2]);
    Amul = contract(Amul, [3 5], Hcen2, [2 3]);
    Amul = contract(Amul, [2 5], Hright, [2 3], [1 4 2 3]);
    Amul = contract(conj(Anew), [1 2 3 4], Amul, [1 2 3 4]);
    Enew = real(getscalarMK(Amul));
    
    end