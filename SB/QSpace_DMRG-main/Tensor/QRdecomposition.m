function [Q,R]= QRdecomposition(T,inputind,varargin)
    totalrankarray = 1:rank(T);
    remainind = totalrankarray(~ismember(totalrankarray, inputind));
    
    [inind, outind] = getlegdir(T);
    [allinT, X] = legflip(T, outind);
    [rank2T, inputId]  = legfuse(allinT, inputind, 'in');
    [rank2T, remainId] = legfuse(rank2T, 1:numel(remainind),'in');
    for itoutind = 1:numel(outind)
        if ismember(outind(itoutind), inputind)
            Xind = find(inputind == outind(itoutind));
            inputId = contract(inputId, Xind, X(itoutind), 2, [1:Xind-1, numel(inputind)+1, Xind:numel(inputind)]);
        elseif ismember(outind(itoutind), remainind)
            Xind = find(remainind == outind(itoutind));
            remainId = contract(remainId, Xind, X(itoutind), 2, [1:Xind-1, numel(remainind)+1, Xind:numel(remainind)]);
        else
            error('strange index!')
        end
    end
    Q = legflip(rank2T,2);
    R = rank2T;
    for i = 1:numel(rank2T.data)
        [q,r] = qr(rank2T.data{i});
        Q.data{i} = q;
        R.data{i} = r;
    end
    Q = contract(conj(inputId), rank(inputId), Q, 1);
    R = contract(R, 2, conj(remainId), rank(remainId));

    if ~isempty(varargin)
        if length(varargin{1})>0 && varargin{1}(end) ~= '*'
            error('please check the varargin itag')
        else
            Q.info.itags{rank(Q)} = varargin{1};
            R.info.itags{1} = strip(varargin{1},'right','*');
        end
    end
