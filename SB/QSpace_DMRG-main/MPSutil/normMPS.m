function a =  normMPS(M,Nsite)
    % <Describe>

    % Calculate norm of MPS M.

    % <Input>

    % M : MPS
    % Nsite

    % <Output>

    % a : [QSpace tensor] 2-rank tensor. Trace of a is norm.

    % Written by Subin Kim (01.02.2025)
    for itN = 1:Nsite
        if itN == 1
            Mval = updateLeft([],[],M(itN),[],[],M(itN));
        elseif itN == Nsite
            Mval = updateLeft(Mval,2,M(itN),[],[],M(itN));
        else
            Mval = updateLeft(Mval,2,M(itN),[],[],M(itN));
        end
    end
    Mval = permute(Mval,[2 1]);
    a = Mval;
    % for i = 1:length(Mval.data)
    %     eig(Mval.data{i})
    % end
end