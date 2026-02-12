function [outputT, X] = legflip(T, ind, varargin)

% < Description >
%
% [output, X] = legflip(T,ind, [,'prime'])
%
% flip the leg direction of the tensor T. It conserves the original leg
% direction and itags.
%
% < Input >
% T : [QSpace object] QSpace tensor
% ind : [numeric array] flipping leg indices.
%
% < Option >
% 'prime' : not remove the prime itags in the 1j symbol. Please refer the QSpace manual draft
%
% < Output >
% outputT : leg-flipped T
% X : 1j symbols used to flip legs
%
% **NB!**
% legflip(legflip(T, ind), ind) can be "-"T !!! due to the QSpace 1j symbol convention
%
% Written by Minsoo Kim (Dec. 19, 2023)

noprime = true;
if ~isempty(varargin)
    if varargin{1} == "prime"
        noprime = false;
    else
        error('please type prime for legflip')
    end
end


X = QSpace(1,numel(ind));
rankT = rank(T);
outputT = T;

for itind=1:numel(ind)
    X(itind) = getIdentity(outputT,ind(itind),'-0');
    if noprime
        X(itind).info.itags = {getitags(X(itind), 1), getitags(X(itind), 1)};
    end
    if ~islegin(T, ind(itind)) %length(itag{ind(itind)})~=0 && itag{ind(itind)}(end)=='*' %out 
        outputT=contract(outputT,ind(itind),X(itind),1, [(1:ind(itind)-1),rankT,(ind(itind):rankT-1)]);
    else %in
        X(itind)=conj(X(itind));
        outputT=contract(outputT,ind(itind),X(itind),1,[(1:ind(itind)-1),rankT,(ind(itind):rankT-1)]);
    end
end


end