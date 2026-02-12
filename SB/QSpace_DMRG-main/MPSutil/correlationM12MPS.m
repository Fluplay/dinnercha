function OO = correlationM12MPS(M1,M2, O, N,varargin)
    % <Description>

    % Overlap <M1|O_i O_i+N|M2> calculated

    % <Input>
    % M1, M2 : MPS (regardless of convention)
    % O : Operator
    % N : chain length calculated

    % <Output>
    % OO : [1x(L-N) numeric] correlation result for each possible lattice point

    % <Option>
    % Selective : select lattice pair calculated. Select left points of calulated pairs.
    %
    % Written by Subin Kim (12.27.2024)

    L = length(M1);
    OO = zeros(1,L-N);
    itNarray = (1:L-N);
    
    while ~isempty(varargin)
        switch varargin{1}
            case 'Selective'
                itNarray = varargin{2};
                OO = zeros(1,length(itNarray));
                varargin(1:2) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end

    index = 1;
    for itN = itNarray
        for i = 1:itN-1
            if i == 1
                T = updateLeft([],[],M1(i),[],[],M2(i));
            else
                T = updateLeft(T,2,M1(i),[],[],M2(i));
            end
        end
        
        if itN ~=1
            T = updateLeft(T,2,M1(itN),O,3,M2(itN));
        else
            T = updateLeft([],[],M1(itN),O,3,M2(itN));
        end

        for i = itN+1:itN+N-1
            T = updateLeft(T,3,M1(i),[],[],M2(i));
        end

        T = updateLeft(T,3,M1(itN+N),Hconj(O) ,3,M2(itN+N));

        for itN2 = ((itN+N+1):L)
            T = updateLeft(T,2,M1(itN2),[],[],M2(itN2));
        end
        OO(index) = trace(T);
        index = index+1;
    end

end