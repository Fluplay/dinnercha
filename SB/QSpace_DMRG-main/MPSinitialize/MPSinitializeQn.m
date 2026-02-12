function Minit = MPSinitializeQn(MPO,S,T, varargin)
    % Function for orbitalHeisenberg
    % < Description >
    %
    %
    % Initialize MPS with D = 1, MPS quantum number is S, T. For model SU(2)xSU(2), charge included or excluded.
    %
    % < Input >
    % MPO : [1 x N QSpace array] MPO of the Hamiltonian.
    %       Each MPO(n) is a rank-4 tensor acting on site n. The order of legs
    %       of MPO(n) is bottom-top-left-right, where the bottom (top) leg
    %       contracts to the physical leg of bra (ket) tensor. 
    % S, T : [numeric] Quantum number of spin/orbital required.
    %
    % < Option >
    % 'itagsoff' : If you use the option 'itagsoff', the itags of the first and second leg of Minit tensors are emptied.
    % 'noCharge' If you exclude charge.
    %
    % < Output >
    % Minit : [1 x N QSpace array] The result  MPS which is obtained.
    %
    % Written by Subin Kim (01.02.2025)
    
    


    itagson = true;
    noCharge = false;
    while ~isempty(varargin)
        switch varargin{1}
            case 'itagsoff'
                itagson = false;
                varargin(1) = [];
            case 'noCharge'
                noCharge = true;
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

    Sarray = generateQarray(Nsite,S);
    Tarray = generateQarray(Nsite,T);
    % Tarray = zeros(Nsite+1,1);

    
    Minit = QSpace(Nsite,1);
    for itN = 1:Nsite
        if itagson
            Anow = getIdentity(Aprev,2,MPO(itN),2,[1 3 2],strcat('bond',char(string(itN))));
        else
            Anow = getIdentity(Aprev,2,MPO(itN),2,[1 3 2]);
        end
        for i = 1:numel(Anow.data)
            Anow.data{i} = Anow.data{i}(:,1,:);
            % Anow.data{i}(1,1) = 1;
        end
        Q = Anow.Q{2}(1,1);
        Snum = Sarray(itN+1);
        Tnum = Tarray(itN+1);
        if noCharge
            Anow = getsub(Anow,[Snum Tnum],2);
        else
            Anow = getsub(Anow,[Q Snum Tnum],2);
        end
        Aprev = Anow;
        % legflip(Aprev,2);
        if itN == Nsite
            Aprev = legflip(Aprev,2);
            if itagson
                Aprev.info.itags{2} = 'end';
            end
        end
        Minit(itN) = Aprev;
    end

    Minit = canonFormMK(Minit,Nsite,[]);
    
    toc2(tobj,'-v');
    
    end