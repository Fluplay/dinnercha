clear

%addpath('/home/dinnercha/local/MuNRG')
%startup
num_threads_SL(8);
U = 2
mu = -U/2
D = 2
ocont = getAcont(0,0,0,0);


% Bethe lattice
[ocont,RhoV2in] = get_Bethe(ocont,D);

T = 1e-8;

Lambda = 2;

N = max(ceil(-2*log(T/500)/log(Lambda)),20);

nz = 2;
Nkeep = 300;

numCorr = 3;

nrgdata = go('data/mytest'); %Work folder

[FF,ZF,SF,IF] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1)

[FF,ZF,SF,EF] = setItag('s00','op',FF(:),ZF,SF(:),IF.E);

NF = quadOp(FF,FF,[]);


HU = (U/2) * sum(NF) * (sum(NF)-EF);

Hmu = sum(NF) * mu;

vac = setItag('L00',getvac(EF)) ;%Vaccum Identity

A0 = getIdentity(vac,2,EF,2,'K00*',[1,3,2]);

HU = contract(A0,'!2*',{HU,'!1',A0});
Hmu = contract(A0,'!2*',{Hmu,'!1',A0});

O_ch = sum(NF/2)

O_sp = SF(end)

[ocont,SE,nrgdataz,ImSEle] = impSE(Hmu,HU,A0,FF,ZF,RhoV2in,T,Lambda,nz,Nkeep,nrgdata,'numCorr',numCorr);


% For Spectral function
SE_H = zeros(1,nz)

ndz = zeros(1,nz);
szz = zeros(1,nz);
for itz = 1:nz
	HU = (U/2) * sum(NF) * (sum(NF)-EF);
	QF = QSpace(size(FF)) %QF = [FF,HU] Kugler IE위해서 ...
	for ito = (1:numel(FF))
		QF(ito) = contract(FF(ito),'!1',HU,[1 3 2]) - contract(HU,'!1',FF(ito));

	end

	QFFF = QSpace(size(FF))

	for ito = (1:numel(FF))
		QFFF(ito) = (contract(QF(ito),'!1',FF(ito),'*')+contract(FF(ito),'!2*',QF(    ito)))/trace(getIdentity(FF(ito),3));
	end



	SE_H(itz) = getEpVal(nrgdataz{itz},QFFF(1))

	QFz = QF(1) - SE_H(itz)*FF(1)
	Op1f = [FF(1);QFz	;FF(1);QFz];
	Op2f = [FF(1);FF(1);QFz  ;QFz];
	zflagf = ones(numel(Op1f),1);
	cflagf = ones(numel(Op1f),1);


	ndz(itz) = getEpVal(nrgdataz{itz},sum(NF));
	Op1c = (sum(NF)-ndz(itz)*EF)/2
	Op2c = Op1c;

	Op1s = SF(end);
	Op2s = Op1s;
	zflagc = zeros(numel(Op1c),1);
	cflagc = -ones(numel(Op1c),1);
	zflags = zeros(numel(Op1s),1);
	cflags = -ones(numel(Op1s),1);
	Ops1 = [Op1f;Op1c;Op1s];
	Ops2 = [Op2f;Op2c;Op2s];
	zflag = [zflagf;zflagc;zflags];
	cflag = [cflagf;cflagc;cflags];
	Adiscs = cell(4+2,nz); % discrete data of two-point correlators (4 fermionic, 2 bosonic)
	[odisc,Adiscs(:,itz),sigmak] = getAdisc(nrgdataz{itz},Ops1,Ops2,ZF,'zflag',zflag,'cflag',cflag);
end



for ita = (1:size(Adiscs,1))
    if ita <= numel(Op1f) % fermionic correlators
        [ocont,Aconts{ita}] = getAcont(odisc,mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3),sigmak,T/5,'alphaz',1/nz);
    else % bosonic correlators
		[ocont,Aconts{ita}] = getAcont(odisc,sum(mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3),2),log(Lambda),T/5,'alphaz',1/nz,'Hfun','CLG');
 	end
end
