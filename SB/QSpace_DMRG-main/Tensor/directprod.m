function AB = directprod(A, B)
% < Description >
%
% AB = directprod(A, B)
%
% direct product two tensors A, B.
%
% < Input >
% A, B : [QSpace object] QSpace tensor
%
% < Option >
%
% < Output >
% AB : [QSpace object rank - rank(A)+rank(B)] direct product of A, B. The leg order is A~, B~.
%
% Written by Minsoo Kim (Dec. 20, 2023)

A3 = appendSingletons(A, '+');
B3 = appendSingletons(B, '-');
AB = contract(A3, rank(A3), B3, rank(A3));

end
