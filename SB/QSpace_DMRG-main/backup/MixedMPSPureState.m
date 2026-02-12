function M = MixedMPSPureState(M)
    % erase uneffective DOF of RMT for each sectors.
    % Used for mixed MPS

    Nsite = numel(M);
    near0 = 1e-40;
    n = numel(M(Nsite).data);
    outLegSpace = cell(n,2);
    colKeepTot = cell(n,1);
    warned = false;
    totEffectiveStates = 0;
    % Fill outLegSpace :  {i,1} is effective outleg number, (i,2) is 2nd leg Quantum number
    for i = 1:n
        V = M(Nsite).data{i};
        colNorm = vecnorm(V);
        colKeep = colNorm > 1e-3;
        colKeepTot{i} = colKeep;
        outLegSpace{i,1} = sum(colKeep(:));
        totEffectiveStates = totEffectiveStates + sum(colKeep(:));
        if totEffectiveStates >= 2
            fprintf('Warning :  Mixed State\n')
            warned = true;
        end

        outLegSpace{i,2} = M(Nsite).Q{2}(i,:);
    end
    % outLegSpace
    % colKeepTot
    for i = 1:n
        V = M(Nsite).data{i};
        colKeep = colKeepTot{i};
        if ~any(colKeep)
            colKeep(1) = 1;
            for j = 1:n
                if isequal(outLegSpace{j,2},M(Nsite).Q{2}(i,:))
                    if outLegSpace{j,1} > sum(colKeep)
                        colKeep(1:outLegSpace{j,1}) = ones(1,outLegSpace{j,1});
                    end
                end
            end
        else

        end
        M(Nsite).data{i} = V(:,colKeep);
        for j = 1:n
            if isequal(outLegSpace{j,2},M(Nsite).Q{2}(i,:))
                if outLegSpace{j,1} > size(M(Nsite).data{i},2) 
                    M(Nsite).data{i}(:,size(M(Nsite).data{i},2)+1 : outLegSpace{j,1}) =...
                        zeros(size(M(Nsite).data{i},1),outLegSpace{j,1}-size(M(Nsite).data{i},2)) + near0;
                end
            end
        end
    end
end