function Mresult = divideOrbitalLeg(M,varargin)
    % <Description>

    % Applied to MPS with physical bond quantum number [U(1)xSU(2)xSU(2)]= [-1 1 1]
    % Divide physical leg to two legs; [-1 1 0] and [0 0 1]
    % Optional : Applied to MPS with physical bond quantum number [SU(2)xSU(2)]= [1 1]
    % Optional : Divide physical leg to two lges; [1 0] and [0 1] 

    % Written by Subin Kim (01.02.2025)

    noCharge = false;

    while ~isempty(varargin)

        switch varargin{1}
            case 'noCharge'
                noCharge = true;
                varargin(1) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end
    % div1, div2: index of quantum number to divide
    if ~noCharge
        Nsite = numel(M);
        Mresult = QSpace(Nsite,1);
        [~,~,Op1,Is]=getLocalSpace('FermionS','Acharge,SU2spin','NC',2);
        Op1 = addSymmetry(Op1,'SU2','q',0);
        Op1.data{1} = 1;
        Op1 = getsub(Op1,[-1 1 0],1);
        [~,~,Op2,Is]=getLocalSpace('FermionS','Acharge,SU2spin','NC',1);
        Op2 = addSymmetry(Op2,'SU2','q',0,'pos',2);
        Ic = getIdentity(Op1,1,Op2,1);
        for itN = 1:Nsite
            Ms = M(itN);
            Ic.info.itags{3} = strcat(Ms.info.itags{3},'*');
            Ic.info.itags{2} = Ms.info.itags{3};
            Ic.info.itags{1} = Ms.info.itags{3};
            if rank(Ms) == 3
                Ms = contract(Ms,3,Ic,3,[1 2 3 4]);
            elseif rank(Ms) == 4
                Ms = contract(Ms,3,Ic,3);
                Ms = contract(Ms,3,Ic,3);
            end
            Mresult(itN) = Ms;
        end
    else
        Nsite = numel(M);
        Mresult = QSpace(Nsite,1);
        [Op1,~]=getLocalSpace('Spin',1/2);
        Op1 = addSymmetry(Op1,'SU2','q',0);
        % Op1 = getsub(Op1,[-1 1 0],1);
        [Op2,~]=getLocalSpace('Spin',1/2);
        Op2 = addSymmetry(Op2,'SU2','q',0,'pos',1);
        Ic = getIdentity(Op1,1,Op2,1);
        for itN = 1:Nsite
            Ms = M(itN);
            Ic.info.itags{3} = strcat(Ms.info.itags{3},'*');
            Ic.info.itags{2} = Ms.info.itags{3};
            Ic.info.itags{1} = Ms.info.itags{3};
            if rank(Ms) == 3
                Ms = contract(Ms,3,Ic,3,[1 2 3 4]);
            elseif rank(Ms) == 4
                Ms = contract(Ms,3,Ic,3);
                Ms = contract(Ms,3,Ic,3);
            end
            Mresult(itN) = Ms;
        end
    end
    end