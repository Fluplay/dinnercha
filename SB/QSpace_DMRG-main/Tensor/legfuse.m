function [fuseT, isometry] = legfuse(T, ind, inorout)

% < Description >
%
% [fuseT, isometry] = legfuse(T, ind, inorout)
%
% fusing tensor legs which have the same direction. The fused leg becomes
% the new last leg.
%
% < Input >
% T : [QSpace object] QSpace tensor
% ind : [numeric array] leg indices want to fuse
%
% < Option >
%
% < Output >
% fuseT : leg fused tensor
% isometry : the isometry used to fuse legs. If one contract fuseT with
% conj(isometry), the T is restored. The last leg is fused leg. 
%
% Written by Minsoo Kim (Dec. 21, 2023)


if inorout == "in"
    legin = true;
elseif inorout == "out"
    legin = false;
else
    error('please type in or out for fuse leg direction')
end

%sanity check
for itind = ind
    if islegin(T, itind) ~= legin
        error('please match the leg direction of indices and inorout')
    end
end

%find the isometry
if numel(ind) == 1
    isometry = getIdentity(T, ind(1));
elseif numel(ind) > 1
    isometry = getIdentity(T, ind(1), T, ind(2));
    for itind = 3:numel(ind)
        id = getIdentity(isometry, itind, T, ind(itind));
        isometry = contract(isometry, itind, id, 1);
    end
end

%fuse the leg
if ~legin
    fuseT = contract(T, ind, isometry, (1:numel(ind)));
else
    isometry = conj(isometry);
    fuseT = contract(T, ind, isometry, (1:numel(ind)));
end


end