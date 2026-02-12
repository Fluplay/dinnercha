function M = MakeSameSpaceMPS(M,M0)
    % Project M into M0 final leg space. Use only if M has bigger domain than M0, such as when M0 is a result of DMRG with initial MPS M. This does not guarantee conservation of physics of M.
    Nsite = numel(M);
    n = numel(M(Nsite).data);
    n0 = numel(M0(Nsite).data);
    for i = 1:n0
        hasSector = false;
        for j = 1:n
            if isequal(M0(Nsite).Q{1}(i,:),M(Nsite).Q{1}(j,:)) && isequal(M0(Nsite).Q{2}(i,:),M(Nsite).Q{2}(j,:)) && isequal(M0(Nsite).Q{3}(i,:),M(Nsite).Q{3}(j,:))
                hasSector = true;
                colKeepNum = size(M0(Nsite).data{j},2);
                M(Nsite)
                M(Nsite).data{i} = M(Nsite).data{i}(:,1:colKeepNum,:);
            end
        end
        if ~hasSector
            frpintf('Warning : uncoverd sector on initial M')
        end
    end
end