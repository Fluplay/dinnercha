function [outputMPO, rank4cell] = getMPO(opmatcell, totalid, varargin)

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

% <Update>
% Updated to deal with various MPOs. totalid can be required to attach 3rd, 4th leg for consistensy.
% Updated by Subin Kim (01.02.2025)
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
    % Fill in rank4cell 
    for itcol = 1:ncol
        for itrow = 1:nrow
            if ~isempty(opmatcell{itrow, itcol}) && isempty(rank4cell{itrow, itcol})
                if rank(opmatcell{itrow, itcol}) >= 3
                    rank4cell{itrow, itcol} = makerank4(opmatcell{itrow, itcol}); %lower, upper, left, right
                else
                    rank4cell{itrow, itcol} =opmatcell{itrow, itcol};
                end
            end
        end
    end


    % finding the information of legs for each row and column to fill the central part.
    % Leave it empty if all elements are zero!
    rowid = QSpace(nrow,1); % left leg
    colid = QSpace(ncol,1); % right leg
    
    for itrow=1:nrow
        rowidFound = false;
        rowEmpty = true;
        for itcol=1:ncol
            if ~isempty(rank4cell{itrow, itcol})
                rowEmpty = false;
                rank2Tensor = rank4cell{itrow,itcol};
            end
            if ~isempty(rank4cell{itrow, itcol}) && rank(rank4cell{itrow, itcol})==4
                rowid(itrow) = getIdentity(totalid, 2, rank4cell{itrow, itcol}, 3, [1 3 2]); % totalid(in), getIdleg(out), left leg of MPO(in)
                rowidFound = true;
                break
            end
        end
        if ~rowidFound && ~rowEmpty
            rowid(itrow) =getIdentity(totalid, 2, makerank4(rank2Tensor), 3, [1 3 2]);
        end
    end
    
    for itcol=1:ncol
        colidFound = false;
        colEmpty = true;
        for itrow=1:nrow
            if ~isempty(rank4cell{itrow, itcol})
                colEmpty = false;
                rank2Tensor = rank4cell{itrow, itcol};
            end
            if ~isempty(rank4cell{itrow, itcol}) && rank(rank4cell{itrow, itcol})==4
                colid(itcol) = conj(getIdentity(totalid,1,rank4cell{itrow, itcol},4,[3 1 2])); % (conj) getIdleg(in), totalid(out), left leg of MPO(out)
                colidFound = true;
                break
            end
        end
        if ~colidFound && ~colEmpty
            colid(itcol) = conj(getIdentity(totalid,1,makerank4(rank2Tensor),4,[3 1 2]));
        end
    end
    
    
    
    
    for itrow=1:nrow
        for itcol=1:ncol
            if ~isempty(rowid(itrow)) && ~isempty(colid(itcol)) % if one of these empty, leave empty
                [rowid(itrow), colid(itcol)] = dimensionForce(rowid(itrow),2,colid(itcol),1);
                if isempty(opmatcell{itrow, itcol})
                    rank4cell{itrow, itcol} = contract(rowid(itrow),2,colid(itcol),1,[1 3 2 4])*near0;
                elseif rank(opmatcell{itrow, itcol}) == 2
                    rank4id = directprod(totalid, cat(getIdentity(rowid(itrow), 3),getIdentity(colid(itcol),2)));
                    rank4cell{itrow, itcol} = contract(rank4id, 2, opmatcell{itrow, itcol},1, [1 4 2 3]);
                end
            else
                rank4cell{itrow, itcol} = QSpace();
            end
        end
    end
    % for i = 1:size(rank4cell,1)
    %     for j = 1:size(rank4cell,2)
    %         i
    %         j
    %         rank4cell{i,j}
    %     end
    % end
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
    for i = 1:size(outputMPO.data)
        outputMPO.data{i}(abs(outputMPO.data{i})<100*near0) = 0;
    
    end
    
    