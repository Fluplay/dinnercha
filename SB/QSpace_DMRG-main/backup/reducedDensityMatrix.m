function Md = reducedDensityMatrix(M)
    M(end) = legflip(M(end),2);
    Nsite = numel(M);
    Md = QSpace(Nsite,1);
    for itN = 1:Nsite
        if(mod(itN,10)==0)
            fprintf('Memory Usage Check: itN = %i\n', itN);
            chkmem();
        end
        if itN==1
            Mcon = Hconj(M(itN));
            j1 = conj(getIdentity(Mcon,1,'-0'));
            Mcon = contract(j1,1,Mcon,1,[2 1 4 3]);
            Ms = contract(M(itN),[1 4],Mcon,[1 4]);
            I2 = getIdentity(Ms,1,Ms,3);
            I2.info.itags{3} = Ms.info.itags{1};
            Ms = contract(Ms,[1 3],I2,[1 2],[3 1 2]);
            Md(itN,1) = Ms;
            clear Ms
        elseif itN==Nsite
            Mcon = Hconj(M(itN));
            j2 = getIdentity(Mcon,2,'-0');
            Mcon = contract(Mcon,2,j2,1,[4 1 3 2]);
            Ms = contract(M(itN),[4 2],Mcon,[4 2]);
            I1 = conj(getIdentity(Ms,1,Ms,3));
            I1.info.itags{3} = Ms.info.itags{1};
            Ms = contract(I1,[1 2],Ms,[1 3]);
            Md(itN,1) = Ms;
            clear Ms
        else
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
            Ms = contract(Ms,[2 4],I2,[1 2],[1 4 2 3]);
            Md(itN,1) = Ms;
            clear Ms
        end
    end
end
 