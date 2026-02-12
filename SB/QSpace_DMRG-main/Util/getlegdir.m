function [inind, outind] = getlegdir(T)

% < Description >
%
% [inind, outind] = getlegdir(T)
% 
% get the direction of each leg
%
% < Input >
% T : [QSpace object] QSpace tensor.
%
% < Option >
% 
% < Output >
% inind  : [numeric array] leg indices which are in.
% outind : [numeric array] leg indices which are out.
% 
%
% Written by Minsoo Kim (Dec. 21, 2023)

inind = [];
outind = [];

for itlegs = 1:rank(T)
    if islegin(T, itlegs)
        inind(end+1) = itlegs;
    else
        outind(end+1) = itlegs;
    end
end


end