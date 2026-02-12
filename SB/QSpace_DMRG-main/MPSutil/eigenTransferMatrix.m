function [Vex,Dex,Vn,Dn] = eigenTransferMatrix(A,B,nKrylov,tol, varargin)

    % <Description>

    % [Vn,Dn] = eigenTransferMatrix(A,B,nKrylov)
    % Calculate the eivenvalue and eigenvector of the transfer matrix using Krylov iteration.
    % Hn = Qn^dagger*T*Qn, where T is the transfer matrix

    % <Input>
    % A : [Nunit x 1 QSpace tensor] Top Part of the transfer matrix
    % B : [Nunit x 1 QSpace tensor] Bottom Part of the transfer matrix, use conj(B)
    % tol: [numeric] Tolerance for the convergence of the Krylov subspace
    % nKrylov: [interger] Krylov subspace dimension

    % <Option>
    % 'x0' : [rank-2 QSpace tensor] Initial guess of the eigenvector. (Default: random)
    % 'exclude' : Exclude the eigenvalue and eigenvector of the transfer matrix
    % 'right' : [logical] If true, the right environment is calculated. (Default: false)

    % <Output>
    % Vex : [1 x nKrylov QSpace tensor] The approximate eigenvector of the transfer matrix with extreme eigenvalue
    % Dex : [numeric] Extreme eigenvalue of the transfer matrix
    % Vn : [1 x nKrylov QSpace tensor] The approximate eigenvector set of the transfer matrix
    % Dn : [numeric] Eigenvalues set of the transfer matrix

    % <Written by Subin Kim (02.04.2025)>

    x0input = false;
    exclude = false;
    isright = false;
    iterNum = 1;
    takeMinEig = false;
    while ~isempty(varargin)
        if numel(varargin) < 2
            error('ERR: Option should be set by a pair of option name and value.');
        end
        switch varargin{1}
            case 'x0'
                x0 = varargin{2};
                x0input = true;
                varargin(1:2) = [];
            case 'exclude'
                excludeEig = varargin{2};
                eigenKet = varargin{3};
                eigenBra = varargin{4};
                exclude = true;
                varargin(1:4) = [];
            case 'right'
                isright = varargin{2};
                varargin(1:2) = [];
            case 'iter'
                iterNum = varargin{2};
                varargin(1:2) = [];
            case 'min'
                takeMinEig = varargin{2};
                varargin(1:2) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end
    
    I = getIdentity(A(1),1);
    if ~x0input
        x0 = getRandDataSB(I);
        x0 = x0/norm(x0);
    end

    if isright
        A = A(end:-1:1);
        for i = 1:numel(A)
            A(i) = permute(A(i),[2 1 3]);
        end
        B = B(end:-1:1);
        for i = 1:numel(B)
            B(i) = permute(B(i),[2 1 3]);
        end
        if ~x0input
            x0 = permute(x0,[2 1]);    
        end
    end
    
    for iter = 1:iterNum
        Hn = zeros(nKrylov,nKrylov);
        Qn = QSpace(nKrylov,1);
        Vn = QSpace(nKrylov,1);

        Nunit = numel(A);
        if numel(B) ~= Nunit
            error('ERR: Input size mismatch.');
        end
        if nKrylov < 1
            error('ERR: Krylov subspace dimension should be larger than 0.');
        end

        beta = norm(x0);
        Qn(1) = x0/beta;
        for itn = (1:nKrylov)
            % apply transfer matrix
            q = Qn(itn);
            q = multiplyTransferMatrix(A,B,q);
            if exclude
                for j = 1:numel(eigenKet)
                    q = q- excludeEig(j)*getscalarMK(contract(Qn(itn),[1 2],eigenBra(j),[1 2]))*eigenKet(j); % T - \lambda |R><L|
                end
            end
            for j  = 1:itn
                % double precision
                hn = zeros(1,2);
                for it2 = 1:2
                    hn(it2) = getscalarMK(contract(q,[1 2],Hconj(Qn(j)),[2 1]));
                    q = q - hn(it2)*Qn(j);
                end
                Hn(j,itn) = sum(hn);
            end
            qNorm = norm(q);
            if qNorm < tol
                break;
            end
            if itn < nKrylov
                Hn(itn+1,itn) = qNorm;
                Qn(itn+1) = q/Hn(itn+1,itn);
            end
        end
        [V,Dn] = eig(Hn);
        Dn = diag(Dn);
        for itn = 1:nKrylov
            Vn(itn) = V(1,itn)*Qn(1);
            for itK = 2:nKrylov
                Vn(itn) = Vn(itn) + V(itK,itn)*Qn(itK);
            end
        end
        if takeMinEig
            [Dex,ind] = min(real(Dn));
            Vex = Vn(ind);
        else
            [Dex,ind] = max(real(Dn));
            Vex = Vn(ind);
        end
        phase = Vex.data{1}(1,1)/abs(Vex.data{1}(1,1));
        Vex = Vex/phase;
        if abs(phase-1) < tol
            Vex = real(Vex);
        end
        x0 = Vex;
    end
end