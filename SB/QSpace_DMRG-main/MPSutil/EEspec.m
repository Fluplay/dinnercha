function Spec = EEspec(Sv)
    % <Describe>

    % Entanglement Spectrum of DMRG result, including selection of adequate bond tensor for spectrum calculation.
    % The spectrum does not include degeneracy from symmetry.

    % <Input>

    % Sv : [(2*Nsweep) x N QSpace array] Sv(n) contains a vector of singular values on
    %     the bond between sites n-1 and n. Sv(1) is empty.

    % <Output>

    % Spec = [N x 1 numeric array] entanglement spectrum sorted.

    % Written by Subin Kim (01.02.2025)

    [m,n] = size(Sv);
    for i = 1:n
        if ~isempty(Sv(1,i))
            break;
        end
    end
    for j = n:-1:1
        if ~isempty(Sv(1,j))
            break;
        end
    end
    Sv = Sv(:,i:j);
    [m,n] = size(Sv);
    for i = 1:m
        if isempty(Sv(i,1))
            i = i-1;
            break;
        end
    end
    Sb = Sv(i,round(n/2)-1);
    Sb = contract(Sb,2,Sb,'2*');
    Spec = EEspecSingle(Sb);
end