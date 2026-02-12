function HconjT = Hconj(T)

% < Description >
%
% HconjT = Hconj(T)
%
% get Hermitian conjugate tensor for rank 2, 3, 4 tensor
%
% < Input >
% T : [QSpace tensor] QSpace tensor. The rank must be 2, 3, 4
%
% < Option >
%
% < Output >
% HconjT : [QSpace tensor] Hermitian conjucation of T.
%           rank 2 -> flip the leg 1, 2.
%           rank 3 -> flip the leg 1, 2. The leg convention of the operator is physical, physical, virtual.
%           rank 4 -> flip the leg [1, 2] and [3, 4] (for the getMPOMK!).
% Written by Minsoo Kim (Dec. 19, 2023)

if rank(T) == 2
    HconjT = permute(T, [2 1], 'conj');
elseif rank(T) == 3
    HconjT = permute(T, [2 1 3], 'conj');
elseif rank(T) == 4
    HconjT = permute(T, [2 1 4 3], 'conj');
else
    error('please check the rank!')
end

end