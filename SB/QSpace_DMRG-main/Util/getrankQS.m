function rankofT = getrankQS(T)
% < Description >
%
% rankofT = getrankQS(T)
%
% find the rank of QSpace object T as counting the number of itags
%
% < Input >
% T : QSpace tensor
%
% < Output >
% rankofT : the rank of T
%
% Written by Minsoo Kim (Dec. 19, 2023)

itag = T.info.itags;
rankofT = length(itag);

end