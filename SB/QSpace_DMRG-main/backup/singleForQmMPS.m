function Minit = singleForQmMPS(Minit)
    Nsite = numel(Minit);
    for i = 1:numel(Minit(Nsite).data)
        Minit(Nsite).data{i} = Minit(Nsite).data{i}(:,1);
    end
end