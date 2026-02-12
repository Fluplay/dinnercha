function [HL, resNorm, resTensor, HLdiff, hl] = gmres_HL_QSpace(Al, C, h, n,nKrylov, tol,varargin)
    % < Description >
    %
    % HL = gmres_HL_QSpace (Al, h, n)
    
    % Solve left environment of h(n) by GMRES methodl, Arnoldi Iteration.
    % b = x*A , b = hl*(1 - P), x_0 = , r_0 = b - A*x_0 
    
    % < Input >
    % Al : [1 x Nunit QSpace array] Left caninocal tensors array
    % h : [1 x Nunit QSpace array] set of nearest neighboor interaction, if h is not a single tensor. In this case, multiHamiltonian = true.
    % n : [numeric] site index of the Hamiltonian

    % <Option>
    % 'right', .. : [logical] If true, the right environment is calculated. (Default: false)
    % 'x0', .. : [rank-2 QSpace tensor] Initial guess of the solution. (Default: random)

    % < Output >
    % HL : [rank-2 QSpace tensor] : Geometric summation of Transfer matrix series, where T_tot = T(n+2)T(n+3)...T(Nunit)T(1)...T(n+1)

    % Written by Subin Kim (02.04.2025)

    isright = false;
    x0input = false;
    while ~isempty(varargin)
        if numel(varargin) < 2
            error('ERR: Option should be set by a pair of option name and value.');
        end
        switch varargin{1}
            case 'right'
                isright = varargin{2};
                varargin(1:2) = [];
            case 'x0'
                x0 = varargin{2};
                x0input = true;
                varargin(1:2) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end

    Nunit = numel(Al);
    if numel(h) ~= Nunit
        error('ERR: Input size mismatch.');
    end

    % if isright
    %     RL = contract(Hconj(C(circularIndex(n+1,Nunit))),2,C(circularIndex(n+1,Nunit)),1);
    % else
    %     RL = contract(Hconj(C(circularIndex(n+2,Nunit))),1,C(circularIndex(n+2,Nunit)),2);
    % end


    if isright
        for i = 1:Nunit
            Al(i) = permute(Al(i),[2 1 3]);
        end
        Al = Al(end:-1:1);
        Al = [Al(2:end);Al(1)];

        for i = 1:Nunit
            h(i) = permute(h(i),[3 4 1 2]);
        end
        h = h(end:-1:1);

        C = C(end:-1:1);
        for i = 1:Nunit
            C(i) = permute(C(i),[1 2]);
        end

        n = Nunit - n + 1;
    end

    RL = contract(Hconj(C(circularIndex(n+1,Nunit))),1,C(circularIndex(n+1,Nunit)),2);
        
    I = getIdentity(Al(circularIndex(n+2,Nunit)),1);
    % I = I/norm(I);   % I is normalized |0>
    Hn = zeros(nKrylov+1,nKrylov);
    Qn = QSpace(nKrylov+1,1);

    hl = updateLeft([],[],Al(n),h(n),4,Al(n));
    hl = contract(Hconj(Al(circularIndex(n+1,Nunit))),[2 3],hl,[1 3]);
    hl = contract(hl,[2 3],Al(circularIndex(n+1,Nunit)),[1 3]);
    b = hl - I*getscalarMK(contract(hl,[1 2],RL,[1 2]));
    if ~x0input
        x0 = getRandDataSB(I);
        x0 = x0*norm(b)/norm(x0);
        % x0 = 0*I;
    end

    Ax0 = x0;
    for i = n+1+1:n+1+Nunit
        Ax0 = updateLeft(Ax0,2,Al(circularIndex(i,Nunit)),[],[],Al(circularIndex(i,Nunit)));
    end
    Ax0 = x0 + I*getscalarMK(contract(x0,[1 2],RL,[1 2])) - Ax0;

    r0 = b - Ax0;
    beta = norm(r0);
    Qn(1) = r0/beta;
    for itn = (1:nKrylov)
        % apply transfer matrix
        q = Qn(itn);
        for i = n+1+1:n+1+Nunit
            q = updateLeft(q,2,Al(circularIndex(i,Nunit)),[],[],Al(circularIndex(i,Nunit)));
        end
        % apply A, including transfer matrix
        q =  Qn(itn) + I*getscalarMK(contract(Qn(itn),[1 2],RL,[1 2])) - q;
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
        Hn(itn+1,itn) = qNorm;
        Qn(itn+1) = q/Hn(itn+1,itn);
    end

    be1 = zeros(nKrylov+1,1);
    be1(1) = beta;

    yn = lsqlin(Hn, be1);
    res =be1-Hn*yn;
    resNorm = norm(res)/norm(b);

    HL = QSpace();
    for itn = (1:nKrylov)
        HL = HL + yn(itn)*Qn(itn);
    end
    HL = x0 + HL;

    % Residue check
    % HLTL = HL;
    % for i = n+1+1:n+1+Nunit
    %     HLTL = updateLeft(HLTL,2,Al(circularIndex(i,Nunit)),[],[],Al(circularIndex(i,Nunit)));
    % end
    % HLdiff = HL-HLTL;
    % resTensor = b - HL - I*getscalarMK(contract(HL,[1 2],RL,[1 2])) + HLTL;

end