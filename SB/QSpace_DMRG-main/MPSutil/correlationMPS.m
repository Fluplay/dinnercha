function OO = correlationMPS(MPS, O, N,varargin)
    % <Description>

    % Overlap <MPS|O_i O_i+N|MPS> calculated

    % <Input>
    % MPS (left-canonical)
    % O : Operator
    % N : chain length calculated

    % <Output>
    % OO : [1x(L-N) numeric] correlation result for each possible lattice point

    % <Option>
    % Selective : select lattice pair calculated. Select left points of calulated pairs.
    %
    % Written by Subin Kim (12.27.2024)

    L = length(MPS);
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
        % for i = 1:itN-1
        %     if i == 1
        %         T = updateLeft([],[],MPS(i),[],[],MPS(i));
        %     else
        %         T = updateLeft(T,2,MPS(i),[],[],MPS(i));
        %     end
        % end
        
        % if itN ~=1
        %     T = updateLeft(T,2,MPS(itN),O,3,MPS(itN));
        % else
            T = updateLeft([],[],MPS(itN),O,3,MPS(itN));
        % end

        for i = itN+1:itN+N-1
            T = updateLeft(T,3,MPS(i),[],[],MPS(i));
        end

        T = updateLeft(T,3,MPS(itN+N),Hconj(O),3,MPS(itN+N));

        for itN2 = ((itN+N+1):L)
            T = updateLeft(T,2,MPS(itN2),[],[],MPS(itN2));
        end
        OO(index) = trace(T);
        index = index+1;
    end

end