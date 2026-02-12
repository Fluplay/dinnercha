function Cleft = updateLeft(Cleft,rankC,B,X,rankX,A)
% < Description >
%
% Cleft = updateLeft(Cleft,rankC,B,X,rankX,A)
%
% Contract the operator Cleft that act on the Hilbert space of the left
% part of the MPS (i.e., left of a given site) with the tensors B, X, and
% A, acting on the given site.
% The details are the same with the same name function in the TN tutorial lecture.
%
% **NB!**
% possible rank of C and X is [rankC rankX],[2 2; 2 3; 2 4; 3 2; 3 3; 3 4; 4 2; 4 3; 4 4]
%
% < Input >
% Cleft : [QSpace tensor] Rank-2 or 3 or 4 QSpace tensor from the left part of the system. If
%       given as empty (i.e., []), then Cleft is considered as the identity
%       tensor of rank 2. If rank-2, its legs are ordered as bottom-top. If rank-3,
%       bottom-top-right. If rank-4, bottom-top-left-right.
% rankC : [integer] Rank of Cleft. If Cleft is given as [], set rankC as 2.
% B, A : [QSpace tensors] Ket tensors, whose legs are ordered as left-right-bottom
%       (local physical). In the contraction, the Hermitian conjugate
%       (i.e., bra form) of B is used, while A is contracted as it is. This
%       convention of inputting B as a ket tensor reduces extra
%       computational cost of taking the Hermitian conjugate of B.
% X : [QSpace tensor] Local operator with rank 2, 3, or 4. If given as [], then X
%       is considered as the identity. If rank-2, its legs are ordered as
%       bottom-top. If rank-3, bottom-top-right. If rank-4, bottom-top-left
%       -right.
% rankX : [integer] Rank of X. If X is given as [], set rankX as 2.
%
% < Output >
% Cleft : [tensor] Contracted tensor. The tensor network diagrams are in
% the TN tutorial material by S. Lee.
%
% Written by Minsoo Kim (Dec.21,2023)


% sanity check
if isempty(Cleft)
    rankC = 2; % regard as the case of the rank-2 identity
end
if isempty(X)
    rankX = 2; % regard as the case of the rank-2 identity
end

if ~ismember([rankC rankX],[2 2; 2 3; 2 4; 3 2; 3 3; 3 4; 4 2; 4 3; 4 4],'rows')
    error('ERR: Invalid ranks of C and X.');
end

B = conj(B); % take complex conjugate to B, without permuting legs

if isempty(Cleft) && isempty(X)
    Cleft = contract(A,[1 3], B,[1 3], [2 1]);
    return
end

if isempty(Cleft)
    if rankX == 2
        T = contract(A,3,X,2);
        Cleft = contract(T, [1 3], B, [1 3], [2 1]);
    elseif rankX == 3
        T = contract(A,3,X,2);
        Cleft = contract(T, [1 3], B, [1 3], [3 1 2]);
    elseif rankX == 4
        T = contract(A,3,X,2);
        Cleft = contract(T, [1 3], B, [1 3] , [4 1 2 3]);
    end
elseif rankC == 2
    if isempty(X)
        T = contract(A,1,Cleft,2);
        Cleft = contract(T,[2 3], B, [3 1], [2 1]);
    elseif rankX == 2
        T = contract(A,1,Cleft,2);
        T = contract(T,2,X,2);
        Cleft = contract(T,[2 3], B, [1 3], [2 1]);
    elseif rankX == 3
        T = contract(A,1,Cleft,2);
        T = contract(T,2,X,2);
        Cleft = contract(T,[2 3], B, [1 3], [3 1 2]);
    elseif rankX == 4
        T = contract(A,1,Cleft,2);
        T = contract(T,2,X,2);
        Cleft = contract(T,[2 3], B, [1 3], [4 1 2 3]);
    end
elseif rankC==3
    if isempty(X)
        T = contract(A,1,Cleft,2);
        Cleft = contract(T,[2 3], B, [3 1], [3 1 2]);
    elseif rankX == 2
        T = contract(A,1,Cleft,2);
        T = contract(T,2,X,2);
        Cleft = contract(T,[2 4], B, [1 3], [3 1 2]);
    elseif rankX == 3
        T = contract(A,1,Cleft,2);
        T = contract(T,[2 4],X,[2 3]);
        Cleft = contract(T,[2 3], B, [1 3], [2 1]);
    elseif rankX == 4
        T = contract(A,1,Cleft,2);
        T = contract(T,[2 4],X,[2 3]);
        Cleft = contract(T,[2 3], B, [1 3], [3 1 2]);
    end
elseif rankC == 4
    if isempty(X)
        T = contract(A,1,Cleft,2);
        Cleft = contract(T,[2 3], B, [3 1], [4 1 2 3]);
    elseif rankX == 2
        T = contract(A,1,Cleft,2);
        T = contract(T,2,X,2);
        Cleft = contract(T,[2 5], B, [1 3], [4 1 2 3]);
    elseif rankX == 3
        T = contract(A,1,Cleft,2);
        T = contract(T,[2 5],X,[2 3]);
        Cleft = contract(T,[2 4], B, [1 3], [3 1 2]);
    elseif rankX == 4
        T = contract(A,1,Cleft,2);
        T = contract(T,[2 5],X,[2 3]);
        Cleft = contract(T,[2 4], B, [1 3], [4 1 2 3]);
    end
end
end