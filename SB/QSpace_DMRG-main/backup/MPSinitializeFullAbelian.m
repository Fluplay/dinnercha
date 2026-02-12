function Minit = MPSinitializeFullAbelian(MPO, varargin)

    % < Description >
    % MPSinitializeFullAbelian which is made to discard minus quantum number
    % [Minit, EGSiter, Eeig, Ulist, Hnowlist, Hmat, finalEK, finalAK] = MPSinitialize(MPO, Nkeep, EorQ, [, 'Enum or Qnum', 'itagsoff'])
    %
    % For given MPO, finding approximate GS MPS by using iterative
    % diagonalization. If Nkeep is larger than the Hilbert space dimension, it
    % becomes exact.
    %
    % < Input >
    % MPO : [1 x N QSpace array] MPO of the Hamiltonian.
    %       Each MPO(n) is a rank-4 tensor acting on site n. The order of legs
    %       of MPO(n) is bottom-top-left-right, where the bottom (top) leg
    %       contracts to the physical leg of bra (ket) tensor. 
    % Nkeep : [numeric] Maximum bond dimension of the MPS to keep.
    % EorQ  : [string] When we find the GS in the last step, we can find the GS
    %       in the total Hilbert space (use 'E') or in the desired quantum number
    %       subspace (use 'Q'). The string 'E' accompanies the number of kept states
    %       for the final step (default=1). The string 'Q' accompanies the
    %       quantum number. There is no default value. The input Qnum must be exist in the varargin.
    %
    % < Option >
    % 'Enum or Qnum' [numeric array] : Refer the description of EorQ
    % 'itagsoff' : If you use the option 'itagsoff', the itags of the first and second leg of Minit tensors are emptied.
    % 'Each' " If you use the option 'Each', resulting MPS is 'Enum'th excited state, excluding other energise.
    %
    % < Output >
    % Minit : [1 x N QSpace array] The result GS MPS which is obtained by
    %       iterative diagonalization. The directions of the first leg of the MPS(1) and the second leg of
    %       the MPS(end) are in and the itags are 'start' and 'end'.
    %       The MPS has left-canonical form.
    %       The leg information of nth MPS is left(bond n-1, in), right(bond n*, out), bottom (s n, in).
    % EGSiter : [numeric] The GS energy.
    % Eeig : [QSpace] Array of eigen values. The number of elements is
    %        the same as Enum (EorQ=='E') or 1 (EorQ=='Q')
    % Ulist : [1 x N QSpace array] Ulist(n) contains eigenstates of each step (lattice site).
    % Hnowlist : [1 x N QSpace array] Hnowlist(n) contains Hamiltonians of each step (lattice site).
    % Hmat : [rank 2 QSpace tensor] Hamiltonian of total system.
    % finalEK : [QSpace] The eigenenergies of the final step. The number of elements is
    %           the same as Enum (EorQ=='E') or 1 (EorQ=='Q'). Empty QSpace when 'Each'.
    % finalAK : [QSpace] The eigenstates of the final step. The dimension of the leg 2 is
    %           the same as Enum (EorQ=='E') or 1 (EorQ=='Q')
    %
    % Written by Minsoo Kim (Dec. 06, 2022)
    % Updated by Minsoo Kim (Dec. 27, 2023): Using updateLeft to update Hprev
    % Updated by Sobun Kim(July. 13, 2024) : add 'Each' argument for EorQ == 'E'
    
    itagson = true;
    
    while ~isempty(varargin)
        switch varargin{1}
            case 'itagsoff'
                itagson = false;
                varargin(1) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end
    
    
    
    
    tobj = tic2;
    
    Nsite = numel(MPO);
    vac = getvac(MPO(1));
    Aprev = vac;
    if itagson
        Aprev.info.itags = {'start','start*'};
    end
        
    Minit = QSpace(Nsite,1);
    for itN = 1:Nsite
        if itagson
            Anow = getIdentity(Aprev,2,MPO(itN),2,[1 3 2],strcat('bond',char(string(itN))));
        else
            Anow = getIdentity(Aprev,2,MPO(itN),2,[1 3 2]);
        end
        for i = 1:numel(Anow.data)
            Anow.data{i} = Anow.data{i}(:,1);
            Qnums = Anow.Q{2}(i,:);
            minusQ = false;
            for n = Qnums
                if n < 0
                    Anow.data{i}(1,1) = 0;
                    minusQ = true;
                    break
                end
            end
            if ~minusQ
                Anow.data{i}(1,1) = 1;
            end
        end
        Ic = getIdentity(Anow,3);
        Anow = contract(Anow,3,Ic,2);
        % legflip(Aprev,2);
        if itN == Nsite
            Aprev = legflip(Aprev,2);
            if itagson
                Aprev.info.itags{2} = 'end';
            end
        end
        Minit(itN) = Aprev;
    end
    
    toc2(tobj,'-v');
    
    end