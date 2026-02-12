function MPO = MajumdarGhoshMPO(J, Nsite)
itagson = true;
while ~isempty(varargin)
    switch varargin{1}
        case 'itagsoff'
            itagson = false;
            varargin(1) = [];
        otherwise
            error('ERR: Unknown input.');
    end
end
[S, info] = getLocalSpace('Spin', 1/2);
I = info.E;

Hloc = cell(4,4);
Hloc{1,1} = I;
Hloc{2,1} = Hconj(S);
Hloc{3,1} = 1e-40*Hconj(S);
Hloc{end,1} = 1e-40*I;
Hloc{end,end} = I;
Hloc{3,2} = I;
Hloc{end,2} = J*S;
Hloc{end,3} = J/2*S;
MPOloc = getMPOMK(Hloc,I);

MPO = QSpace(Nsite,1);
MPO(1) = getMPOMK(Hloc,I, "start");
MPO(end) = getMPOMK(Hloc,I, "end");

for itN = 2:Nsite-1
    MPO(itN) = MPOloc;
end

MPO = attachitags2MPO(MPO,'s','O');


end