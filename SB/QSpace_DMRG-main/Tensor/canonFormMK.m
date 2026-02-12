function [M,S,dw] = canonFormMK (M,id,Nkeep)

% < Description >
%
% [M,S,dw] = canonFormMK (M,id,Nkeep)
%
% Obtain the canonical forms of MPS, depending on the index id of the
% target bond. The left part of the MPS, M(1), ..., M(id), is brought into
% the left-canonical form and the right part of the MPS; M(id+1), ...,
% M(end), into the right-canonical form. Thus, if id is 0, the result is
% purely right-canonical form; if id is numel(M), the result is purely
% left-canonical form.  
%       **NB!**
%       All tensor of input MPS M must have the arrow direction in, out, in,
%       except the last tensor. The arrow directdion of the last tensor
%       must be in, in, in.
%
% < Input >
% M : [QSpace array] MPS of length numel(M). Each cell element is a rank-3
%       tensor, where the first, second, and third dimensions are
%       associated with left, right, and bottom (i.e., physical) legs,
%       respectively.
%       **NB!**
%       The direction of the first leg of M(1) and the second leg of M(end) must be in!
% id : [integer] Index for the bond connecting the tensors M(id) and
%       M(id+1). With respect to the bond, the tensors to the left
%       (right) are brought into the left-(right-)canonical form.
% Nkeep : [number] Maximal number of singular values to keep at each SVD.
%       If set empty ([]), it is interpreted as Inf, meaning no truncation
%       by the number of singular values.
%
% < Output >
% M : [QSpace array] Left-, right-, or bond-canonical form from input M,
%       depending on id, as follows:
%       * id == 0: right-canonical form
%       * id == numel(M): left-canonical form
%       * otherwise: bond-canonical form
%       The leg direction is
%       * out in in : right-canonical part
%       * in out in : left-canonical part
% S : [QSpace column vector] Singular values at the bond between M(id) and M(id+1)
%       if 0 < id < numel(M); the norm of the MPS if id = 1 or numel(M).
%       The leg direction is in in
% dw : [column vector] Vector of length numel(M)-1. dw(n) means the
%       discarded weight (i.e., the sum of the square of the singular  
%       values that are discarded) at the bond between M(n) and M(n+1).
%
% Written by Minsoo Kim (Dec. 21, 2023)


% % check the integrity of input
if (numel(id) ~= 1) || (round(id) ~= id)
    error('ERR: 2nd input ''id'' needs to be a single integer.');
elseif (id < 0) || (id > numel(M))
    error('ERR: the 2nd input ''id'' needs to be in a range (0:numel(M))');
end

if ~islegin(M(1), 1)
    error('ERR: diection of the first leg of MPS(1) must be in');
elseif ~islegin(M(end), 2)
    error('ERR: diection of the second leg of MPS(end) must be in');
end

dw = zeros(numel(M)-1,1); % discarded weights

% % Bring the left part of MPS into the left-canonical form
for it = (1:id-1)
    itag = getitags(M(it), 2);
    if ~isempty(itag) && itag(end) == '*'
        itag = itag(1:end-1);
    end
    [M(it),S,Vd, dw(it)] = svdMK(M(it), [1 3], Nkeep, itag);
    M(it) = permute(M(it), [1 3 2]);
    M(it+1) = contract(S*Vd, 2, M(it+1), 1);
end

% % Bring the right part into the right-canonical form
for it = (numel(M):-1:id+2)
    itag = getitags(M(it), 1);
    if ~isempty(itag) && itag(end) == '*'
        itag = itag(1:end-1);
    end
    [U, S, M(it), dw(it-1)] = svdMK(M(it), 1, Nkeep, itag);
    M(it-1) = contract(M(it-1), 2, U*S, 1, [1 3 2]);
end

if id == 0 % purely right-canonical form
    itag = getitags(M(1), 1);
    if ~isempty(itag) && itag(end) == '*'
        itag = itag(1:end-1);
    end
    [U, S, M(1)] = svdMK(M(1), 1, [], itag); % no trucation
    M(1) = legflip(M(1), 1);
elseif id == numel(M) % purely left-canonical form
    itag = getitags(M(end), 2);
    if ~isempty(itag) && itag(end) == '*'
        itag = itag(1:end-1);
    end
    [M(end), S, Vd] = svdMK(M(end), [1 3], [], itag); % no trucation
    M(end) = permute(M(end), [1 3 2]);
    M(end) = legflip(M(end), 2);
else % bond-canonical form
    itag = getitags(M(id+1), 1);
    if ~isempty(itag) && itag(end) == '*'
        itag = itag(1:end-1);
    end
    T = contract(M(id),2,M(id+1),1);
    [M(id),S,M(id+1), dw(id)] = svdMK(T,[1 2],Nkeep, itag);
    M(id) = permute(M(id),[1 3 2]);
end

end