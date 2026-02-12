function inorout = islegin(T, legind)

% < Description >
%
% inorout = islegin(T, legind)
% 
% find the direction of the leg. in: true, out: false
%
% < Input >
% T : [QSpace object] QSpace tensor.
% legind : [numeric] leg index which we want to know the direction.
%
% < Option >
% 
% < Output >
% inorout : [Bullean] leg direction. in: true, out: false
% 
%
% Written by Minsoo Kim (Dec. 21, 2023)

itag = getitags(T, legind);
if isempty(itag) || itag(end) ~= '*'
    inorout = true;
elseif itag(end) == '*'
    inorout = false;
else
    error('can not determine the leg direction')
end


end