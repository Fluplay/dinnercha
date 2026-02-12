function [Anew, Enew, resNorm] = eigs_HC_QSpace (Al,h,Ar,HL,HR,Aold,nKrylov,tol)
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
        Amul1 = updateLeft([],[],Al,h,4,Al);
        Amul1 = contract(Amul1,2,As(itn),1,[1 4 2 3]);
        Amul1 = contract(Amul1,[2 4],Ar,[1 3]);
        Amul1 = contract(Amul1,[2 3],conj(Ar),[3 2]);
        Amul2 = contract(HL,2,As(itn),1);
        Amul3 = contract(As(itn),2,HR,2);
        Amul = Amul1+Amul2+Amul3;
        alphas(itn) = real(getscalarMK(contract(conj(As(itn)),(1:2),Amul,(1:2))));

        % if itn == 1
        %     fprintf('initial eig is %f\n',alphas(1));
        % end

        cnt = cnt+1;
        if itn < nKrylov
            for it2 = (1:2) % do twice to further reduce numerical noise
                for itK = 1:itn
                    T = contract(conj(As(itK)), [1 2] ,Amul, [1 2]);
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

    end
    
    Hkrylov = diag(betas(1:cnt-1),-1);
    Hkrylov = Hkrylov + Hkrylov' + diag(alphas(1:cnt));
    
    [V,D] = eig(Hkrylov);
    [Enew,minid] = min(diag(D));
    Anew=V(1,minid)*As(1);
    for itK = 2:cnt
        Anew=Anew+V(itK,minid)*As(itK);
    end
    % Enew

    Amul1 = updateLeft([],[],Al,h,4,Al);
    Amul1 = contract(Amul1,2,Anew,1,[1 4 2 3]);
    Amul1 = contract(Amul1,[2 4],Ar,[1 3]);
    Amul1 = contract(Amul1,[2 3],conj(Ar),[3 2]);
    Amul2 = contract(HL,2,Anew,1);
    Amul3 = contract(Anew,2,HR,2);
    Amul = Amul1+Amul2+Amul3;
    Enew = real(getscalarMK(contract(Amul,[1 2],conj(Anew),[1 2])));

    resNorm = norm(Enew*Anew - Amul)/abs(Enew);

    
    end