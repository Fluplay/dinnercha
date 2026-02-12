function M = bondDimDown(M,D,div)
    % <Description>
    
    % Decreasing bond dimension of certain MPS. It is recommended to use DMRG insted of this 
    % when high-quality MPS is required.
    %
    % D : maximum bond Dimension
    % div : the number of iteration of gradually decreasing bond dimension
    % Writted by Subin Kim (12.27.2024)
    Nsite = numel(M);
    Dmax = 0;
    for i = 1:Nsite
        dimM = dim(M(i));
        Dmax = max(Dmax,dimM(2));
    end
    if Dmax <= D
        return;
    end
    Darray = round(linspace(Dmax,D,div));
    for i = 2:length(Darray)
        M = canonFormSB(M,Nsite,Darray(i));
    end
end