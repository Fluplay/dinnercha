function [Al, Ar, C, err, Almin,Armin,Cmin, errMin, errHist] = VUMPS_GS_QSpace (Al, Ar,C,h,errTol,varargin)
    % VUMPS with nearest interaction case, sequential multicell update, without dynamic bond control
    % Case without Controlled bond expansion bond control, MPO Hamiltonian is not updated yet.
    % % Multicell is implied in serial method.

    % < Input >
    % Al, Ar, C : [1 x Nunit QSpace array] Initial array of unit cell elements
    % h : [rank-4 QSpace Tensor] nearest neighboor interaction
        % 2           4
        % |           |
        % |-----------|
        % |           |
        % 1           3
    % (Optional) h : [1 x Nunit QSpace array] set of nearest neighboor interaction, if h is not a single tensor. In this case, multiHamiltonian = true.
    % errTol : [numeric] error tolerance for convergence

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
    % < Output >
    % Al : [1 x Nunit QSpace array] The multicell left unitary tensor set for variational minimum uMPS.
    % Ar : [1 x Nunit QSpace array] The multicell right unitary tensor set for variational minimum uMPS.
    % Ac : [1 x Nunit QSpace array] The multicell center site tensor set for variational minimum.
    %         Ac(k) = Al(k)C = C(k-1)
    % C : [1 x Nunit QSpace array] The multicell bond tensor set for variational minimum.
    % Eiter : [Nsweep numeric array] Each step energy.
    
    % Written by Subin Kim (Nov, 12,2024)
    tobj = tic2;
    % default vales of optional input parameters
    nKrylovGMRES = 5;
    nKrylovEig = 5;
    tol = 1e-10;
    iterMax = inf;
    iterGMRES = 5;
    iterEig = 5;
    % parse options
    while ~isempty(varargin)
        if numel(varargin) < 2
            error('ERR: Option should be set by a pair of option name and value.');
        end
        switch varargin{1}
            case 'KrylovGMRES'
                nKrylovGMRES = varargin{2};
                varargin(1:2) = [];
            case 'KrylovEig'
                nKrylovEig = varargin{2};
                varargin(1:2) = [];
            case 'iterGMRES'
                iterGMRES = varargin{2};
                varargin(1:2) = [];
            case 'iterEig'
                iterEig = varargin{2};
                varargin(1:2) = [];
            case 'tol'
                tol = varargin{2};
                varargin(1:2) = [];
            case 'iterMax'
                iterMax = varargin{2};
                varargin(1:2) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end
    Nunit = numel(Al);
    if(length(h) == 1)
        repmat(h,Nunit,1);
    elseif (length(h) == Nunit)
    else
        error('ERR: The number of Hamiltonian should be 1 or Nunit.');
    end

    err = inf(Nunit,1);
    errMin = inf;
    iterCount = 0;
    errHist = struct();     % errHist(iterCount, Nunit).{resNormHL,resNormHR, resNormAC,resNormCl, resNormCr, err}
    % resNormHL(Nnit,iterGMRES), resNormAC(1, iterEig), err(1,2) = [errL, errR]
    while max(err) > errTol && iterCount < iterMax
        iterCount = iterCount + 1;
        for itN = 1:Nunit
            % Optional: Dynamic bond control


            % Calculate Explit terms of effective Hamiltonians
            [HLout, HLin, HRin, HRout,errHist(iterCount,itN).resNormHL,errHist(iterCount,itN).resNormHR] ...
            = HeffTerms(Al, Ar, C, h, itN,nKrylovGMRES,tol,iterGMRES); % HLin = HLout*TL, HRin = TR*HLout
            HLout = 0.5*(HLout + Hconj(HLout));
            HLin = 0.5*(HLin + Hconj(HLin));
            HRin =  0.5*(HRin + Hconj(HRin));
            HRout = 0.5*(HRout + Hconj(HRout)); 


            % eigensolvers
            [ACnew,Cl,Cr,errHist(iterCount,itN).resNormAC,errHist(iterCount,itN).resNormCl,errHist(iterCount,itN).resNormCr] = eigsSolver_QSpace(Al,Ar,C,h, HLout,HLin,HRin,HRout,itN,nKrylovEig,tol,iterEig);

            % Polar decomposition
            [Al(itN), errL] = polarDecompositionAl(ACnew,Cr);
            [Ar(itN), errR] = polarDecompositionAr(ACnew,Cl);

            % return tensors
            C(circularIndex(itN-1,Nunit)) = Cl;
            C(itN) = Cr;
            err(itN) = max(errL,errR);
            errHist(iterCount,itN).err = [errL,errR];
            disptime(['Sweep #',sprintf('%i/%i, n = %i',[iterCount,iterMax,itN]),', err = ', ...
                sprintf('(%d, %d)',[errL, errR])]);
        end

        if max(err) < errMin
            errMin = max(err);
            Almin = Al;
            Armin = Ar;
            Cmin = C;
        end

    end
    toc2(tobj,'-v');


    

end
    