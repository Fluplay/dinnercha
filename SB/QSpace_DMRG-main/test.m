
addpath(genpath(pwd))


contract(m(1),'!2',conj(m(1)),'!2')

m = MPSinitializeQnSingle(mpo,10)

L1 = contract(m(1),'!2',conj(m(1)),'!2')

for i = 2:numel(m)-1

    Ls1 = contract(L1,m(i));

    L1 = contract(Ls1,'!2',conj(m(i)),'!2');

end

    

LL = contract(L1,m(end))