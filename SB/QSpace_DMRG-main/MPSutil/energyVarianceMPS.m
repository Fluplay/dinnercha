function  var= energyVarianceMPS(M,MPO,Nsite)
    % <Describe>

    % Calculate energy variance of MPS

    % <Input>

    % M : [1 x Nsite QSpace array] MPS calculated
    % MPO : [1 x Nsite QSpace array] MPO
    % Nsite 

    % <Output>
    % var : [numeric] variance of Hamiltonian

    % Written by Subin Kim (01.02.2025)
    
    H2 = QSpace(1,Nsite);
    for i = 1:Nsite
        O = contract(MPO(i),2,MPO(i),1);
        I1 = getIdentity(O,2,O,5);
        itag = getitags(O,2);
        O = contract(conj(I1),[1 2],O,[2 5]);
        O = setitags(O,1,itag);
        I2 = getIdentity(O,3,O,5);
        itag = getitags(O,3);
        O = contract(O,[3 5],I2,[1 2]);
        O = setitags(O,4,itag);
        O = permute(O,[2 3 1 4]);
        H2(i) = O;
    end
    squareExp = HvalMPS(M,H2,Nsite);
    Hval = HvalMPS(M,MPO,Nsite);
    expSquare = contract(Hval,2,Hval,1);
    var = trace(squareExp)-trace(expSquare);
end