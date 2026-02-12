function [outputMPO, rank4cell] = getMPOMK(opmatcell, totalid, varargin)

% < Description >
%
% [outputMPO, rank4cell] = getMPOMK(opmatcell, totalid, [, 'start or end'])
%
% Making rank 4 MPO from cell of operators. leg direction is lower, upper, left, right.
%
% < Input >
% opmatcell : [M x N cell array] M x N cell array containing operators.
%             Each element are QSpace object. Rank of each tensor is 2~4.
%             You can leave the cell elements empty if there is no operator.
%             **NB!**
%             opmatcell{1,1} and opmatcell{end,end} shoud be the same as totalid. 
%             The rank of opmatcell{end,1} must be 2.
%             The leg direction of Rank 3 tensors in the first col must be in.
%             The leg direction of Rank 3 tensors in the last row must be out.
%
% totalid : [QSpace tensor] The identity of local Hilbert space.
%
% < Option >
% 'start or end' : 'start' -> for the left end
%                  'end' -> for the right end
%
% < Output >
% outputMPO : [QSpace tensor] rank 4 MPO. leg direction is lower, upper, left, right.
% rank4cell : [M x N cell array] M x N cell array containing rank 4 operators.
%
% Written by Minsoo Kim (Dec. 21, 2023)
if ~isempty(varargin)
    if numel(varargin) > 2
        error("varargin: none, start, end");
    end
    if varargin{1} ~= "start" && varargin{1} ~= "end"
        error("varargin: none, start, end");
    end
end

[nrow, ncol] = size(opmatcell);
rank4cell = cell(nrow,ncol);


near0 = 1e-40; % You can change the value to very small value (e.g. 1e-40) if there is a bug!

% add singletons to make rank 4 for not empty elements in the first and last col and row.
% I suggest filling the first column and the last row!
for itedgecol = [1 ncol]
    for itrow = 1:nrow
        if ~isempty(opmatcell{itrow, itedgecol}) && isempty(rank4cell{itrow, itedgecol})
            rank4cell{itrow, itedgecol} = makerank4(opmatcell{itrow, itedgecol}); %lower, upper, left, right
        end
    end
end
for itedgerow = [1 nrow]
    for itcol = 1:ncol
        if ~isempty(opmatcell{itedgerow, itcol}) && isempty(rank4cell{itedgerow, itcol})
            rank4cell{itedgerow, itcol} = makerank4(opmatcell{itedgerow, itcol}); %lower, upper, left, right
        end
    end
end
% fill the blank edge
for itrow = 1:nrow
    if isempty(rank4cell{itrow, 1})
        if ~isempty(rank4cell{end, itrow})
            rank4cell{itrow, 1} = Hconj(rank4cell{end, itrow})*near0;
        elseif ~isempty(rank4cell{itrow, end})
            rank4cell{itrow, 1} = rank4cell{itrow, end}*near0;
        else
            error('please check the opmatcell!')
        end
    end
    if isempty(rank4cell{itrow, end})
        if ~isempty(rank4cell{itrow, 1})
            rank4cell{itrow, end} = rank4cell{itrow, 1}*near0;
        elseif ~isempty(rank4cell{1, itrow})
            rank4cell{itrow, end} = Hconj(rank4cell{1, itrow})*near0;
        else
            error('please check the opmatcell!')
        end
    end
end
for itcol = 1:ncol
    if isempty(rank4cell{end, itcol})
        if ~isempty(rank4cell{itcol, 1})
            rank4cell{end, itcol} = Hconj(rank4cell{itcol, 1})*near0;
        elseif ~isempty(rank4cell{1, itcol})
            rank4cell{end, itcol} = rank4cell{1, itcol}*near0;
        else
            error('please check the opmatcell!')
        end
    end
    if isempty(rank4cell{1, itcol})
        if ~isempty(rank4cell{end, itcol})
            rank4cell{1, itcol} = rank4cell{end, itcol}*near0;
        elseif ~isempty(rank4cell{itcol, end})
            rank4cell{1, itcol} = Hconj(rank4cell{itcol, end})*near0;
        else
            error('please check the opmatcell!')
        end
    end
end

% finding the information of legs for each row and column to fill the central part.
% Leave it empty if all elements are zero!
rowid = QSpace(nrow,1); % left leg
colid = QSpace(ncol,1); % right leg

for itrow=1:nrow
    for itcol=1:ncol
        if ~isempty(rank4cell{itrow, itcol}) && rank(rank4cell{itrow, itcol})==4
            rowid(itrow) = getIdentity(totalid, 2, rank4cell{itrow, itcol}, 3, [1 3 2]); % totalid(in), getIdleg(out), left leg of MPO(in)
            break
        end
    end
end

for itcol=1:ncol
    for itrow=1:nrow
        if ~isempty(rank4cell{itrow, itcol}) && rank(rank4cell{itrow, itcol})==4
            colid(itcol) = conj(getIdentity(totalid,1,rank4cell{itrow, itcol},4,[3 1 2])); % (conj) getIdleg(in), totalid(out), left leg of MPO(out)
            break
        end
    end
end




for itrow=2:nrow-1
    for itcol=2:ncol-1
        if isempty(opmatcell{itrow, itcol})
            rank4cell{itrow, itcol} = contract(rowid(itrow),2,colid(itcol),1,[1 3 2 4])*near0;
        elseif rank(opmatcell{itrow, itcol}) == 2
            rank4id = directprod(totalid, getIdentity(rank4cell{itrow, 1}, 3));
            rank4cell{itrow, itcol} = contract(rank4id, 2, opmatcell{itrow, itcol},1, [1 4 2 3]);
        end
    end
end


% Concatenate all rank4 cells
if ~isempty(varargin) && varargin{1} == "start"
    outputMPO = cat(rank4cell{end,:},4);
elseif ~isempty(varargin) && varargin{1} == "end"
    outputMPO = cat(rank4cell{:,1},3);
else
    colsum = cell(nrow,1);
    for itrow=1:nrow
        colsum{itrow} = cat(rank4cell{itrow,:}, 4);
    end
    outputMPO = cat(colsum{:}, 3);
end

end

