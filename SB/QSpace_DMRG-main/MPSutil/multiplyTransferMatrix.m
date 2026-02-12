function v0 = multiplyTransferMatrix(A,B,v0,varargin)
    isright = false;
    while ~isempty(varargin)
        if numel(varargin) < 2
            error('ERR: Option should be set by a pair of option name and value.');
        end
        switch varargin{1}
            case 'right'
                isright = varargin{2};
                varargin(1:2) = [];
            otherwise
                error('ERR: Unknown input.');
        end
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
    end

    Nunit = numel(A);
    for i = 1:Nunit
        v0 = updateLeft(v0,rank(v0),B(i),[],[],A(i));
    end
end