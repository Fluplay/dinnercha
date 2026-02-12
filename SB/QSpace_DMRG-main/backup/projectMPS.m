function Minit = projectMPS(Minit,Q)
    % project MPS fial bond
    Nsite = numel(Minit);
    Minit(Nsite) = getsub(Minit(Nsite),Q,2);
end