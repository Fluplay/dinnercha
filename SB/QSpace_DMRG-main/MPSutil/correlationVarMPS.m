function OO = correlationVarMPS(MPS, O, N,varargin)
    % <Description>

    % Variation of <MPS|O_i-N O_i|MPS> calculated

    % <Input>
    % MPS (left-canonical)
    % O : Operator
    % N : chain length calculated

    % <Output>
    % OO : [1x(L-N) numeric] variation of correlation

    % <Option>
    % Selective : select lattice pair calculated. Select right points of calulated pairs. (Sorry for convention different from 'correlationM1M2MPS.m')
    %
    % Written by Subin Kim (12.27.2024)
    L = length(MPS);
    OO = zeros(1,L-N);
    itNarray = N + (1:L-N);
    
    C = getIdentity(O,3,O,3);
    Osq = contract(O,2,O,1);
    O2 = contract(Osq,'24',C,'12');
    O2d =Hconj(contract(Osq,'24',C,'21'));

    while ~isempty(varargin)
        switch varargin{1}
            case 'Selective'
                itNarray = N + varargin{2};
                OO = zeros(1,length(itNarray));
                varargin(1:2) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end

    index = 1;
    for itN = itNarray
        T = updateLeft([],[],MPS(itN-N),O,3,MPS(itN-N));
        T2 = updateLeft([],[],MPS(itN-N),O2,3,MPS(itN-N));
        for i = itN-N+1:itN-1
            T = updateLeft(T,3,MPS(i),[],[],MPS(i));
            T2 = updateLeft(T2,3,MPS(i),[],[],MPS(i));
        end
        T = updateLeft(T,3,MPS(itN),Hconj(O) ,3,MPS(itN));
        T2 = updateLeft(T2,3,MPS(itN),O2d ,3,MPS(itN));
        for itN2 = ((itN+1):L)
            T = updateLeft(T,2,MPS(itN2),[],[],MPS(itN2));
            T2 = updateLeft(T2,2,MPS(itN2),[],[],MPS(itN2));
        end
        OO(index) = trace(T2)-trace(T)^2;
        index = index+1;
    end

end