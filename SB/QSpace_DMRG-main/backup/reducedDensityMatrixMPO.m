function Md = reducedDensityMatrixMPO(M)
    M(end) = legflip(M(end),2);
    Nsite = numel(M);
    Md = QSpace(Nsite,1);
    for itN = 1:Nsite
            Mcon = Hconj(M(itN));
            j1 = conj(getIdentity(Mcon,1,'-0'));
            j2 = getIdentity(Mcon,2,'-0');
            Mcon = contract(j1,1,Mcon,1);
            Mcon = contract(Mcon,2,j2,1,[4 1 3 2]);
            Ms = contract(M(itN),4,Mcon,4);
            I1 = conj(getIdentity(Ms,1,Ms,4));
            I1.info.itags{3} = Ms.info.itags{1};
            Ms = contract(I1,[1 2],Ms,[1 4]);
            I2 = getIdentity(Ms,2,Ms,4);
            I2.info.itags{3} = Ms.info.itags{2};
            Ms = contract(Ms,[2 4],I2,[1 2],[2 3 1 4]);
            Ms.info.itags{2} = strcat(Ms.info.itags{1},'*');
            Md(itN,1) = Ms;
            clear Ms
    end
end
 