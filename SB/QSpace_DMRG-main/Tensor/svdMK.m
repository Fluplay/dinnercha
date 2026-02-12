function [U, S, Vd, dw] = svdMK(T, inputind, Nkeep, varargin)

% < Description >
%
% [U,S,Vd,dw] = svdMK(T, inputind, Nkeep [,itags])
%
% Singular value decomposition of tensor such that T = U*diag(S)*Vd. (Note
% that it is not U*S*V' as in the MATLAB built-in function 'svd'.) If
% truncation criterion Nkeep is given, the tensors are truncated
% to discard the smallest singular values and the corresponding singular
% vectors.
%
% < Input >
% T : [tensor] Tensor.
% inputind : [integer vector] Indices of T to be associated with U. For example,
%       if rankT == 4 and idU == [1 3], the result U is rank-3 tensor whose
%       1st and 2nd legs correspond to the 1st and 3rd legs of T. The 3rd
%       leg of U is associated with the 1st leg of diag(S). And Vd is
%       rank-3 tensor whose 2nd and 3rd legs correspond to the 2nd and 4th
%       legs of T. Its 1st leg is associated with the 2nd leg of diag(S).
% Nkeep : [number] Maximal number of singular values to keep. If set empty
%       ([]), it is interpreted as Inf, meaning no truncation by the number
%       of singular values.
%
%< Option >
% itags : The itags of the last leg of U, all leg of S, and the first leg of Vd. It must not have '*'
%
%
% < Output >
% U : [QSpace tensor] Tensor describing the left singular vectors. Its last leg
%       contracts with diag(S). The earlier legs are specified by input
%       idU; their order is determined by the ordering of idU.
% S : [rank-2 QSpace tensor] The column vector of singular values. If there were no
%       truncation, norm(S) indicates the norm of the tensor.
% Vd : [QSpace tensor] Tensor describing the right singular vectors. Its 1st leg
%       contracts with diag(S). The later legs conserve the order of the
%       legs of input T.
% dw : [numeric] Discarded weight (= sum of the square of the singular
%       values truncated).
%
% Written by Minsoo Kim (Dec. 21, 2023)

totalrankarray = 1:rank(T);
remainind = totalrankarray(~ismember(totalrankarray, inputind));

[inind, outind] = getlegdir(T);
[allinT, X] = legflip(T, outind);

[rank2T, inputId]  = legfuse(allinT, inputind, 'in');
[rank2T, remainId] = legfuse(rank2T, 1:numel(remainind), 'in');

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

if isempty(Nkeep) 
    [U, S, Vd, svdinfo] = svdQS(rank2T, 2);
elseif isint(Nkeep)
    [U, S, Vd, svdinfo] = svdQS(rank2T, 2,'Nkeep',Nkeep);
else
    error('Nkeep must be an integer')
end

dw = svdinfo.svd2tr;

U = QSpace(U); S = diag(QSpace(S)); Vd = QSpace(Vd);

U = contract(conj(inputId), rank(inputId), U, 1);
Vd = contract(Vd, 2, conj(remainId), rank(remainId));

if ~isempty(varargin)
    if length(varargin{1})>0 && varargin{1}(end) == '*'
        error('please check the varargin itag')
    else
        U.info.itags{rank(U)} = strcat(varargin{1},'*');
        S.info.itags = {varargin{1}, varargin{1}};
        Vd.info.itags{1} = strcat(varargin{1},'*');
    end
end

%% checking U*S*Vd = T

% restoreT = contract(U, rank(U), S, 1);
% restoreT = contract(restoreT, rank(restoreT), Vd, 1);
% [~, sortind] = sort([inputind, remainind]);
% restoreT = permute(restoreT, sortind);
% diff = restoreT - T;
% svderr = contract(diff, (1:rank(T)), conj(diff), (1:rank(T)))


end