function Spec = EEspecSingle(Sb)
    % <Description>

    % Entanglement Spectrum from single 'reduced density matrix' tensor.
    % Degeneracy from symmetry is excluded.

    % Written by Subin Kim (01.02.2025)
    Sb  = diag(Sb);
    Spec = [];
    for i = 1:numel(Sb.data)
        Spec = [Spec; abs(Sb.data{i})];
    end
    Spec = sort(Spec,'descend');
end