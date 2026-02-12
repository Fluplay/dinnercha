function Af = legfuseRank3(A,varargin)
    % Af: rank-2, [1 2] fuse
    divide = false;
    originForm = QSpace();
    while ~isempty(varargin)
        switch varargin{1}
            case 'divide'
                divide = true;
                originForm = varargin{2};
                varargin(1:2) = [];
            otherwise
                error('ERR: Unknown input.');
        end
    end
    if ~divide
        Af = getIdentity(A,3);
        for i = 1:numel(A.data)
            sizeA = size(A.data{i});
            Af.data{i} = reshape(A.data{i},[sizeA(1)*sizeA(2),sizeA(3:end)]);
        end
    else
        Af = originForm;
        for i = 1:numel(Af.data)
            sizeA = size(A.data{i});
            Af.data{i} = reshape(A.data{i},[sqrt(sizeA(1)),sqrt(sizeA(1)),sizeA(2:end)]);
        end
    end
end
