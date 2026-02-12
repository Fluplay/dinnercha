function get_Bethe_DMFT_2orb(U,mu,iter,numCorr,D,resname,RhoV2init,T,Lambda,Nkeep,nz)
	%%%todo
		% => it should be merged with get_Bethe_DMFT.m
	%%%


	RhoV2_arr = []
	RhoV2in_arr = RhoV2init;
	G_arr = [];
	A_arr = [];
	SE_arr = [];
	Delta_arr = [];
	nd_arr = [];
	dd_arr = [];



	sym = 'Acharge(:),SU2spin';
	[FF,ZF,SF,IF] = getLocalSpace('FermionS',sym,'NC',2);
	[FF,ZF,SF,EF] = setItag('s00','op',FF,ZF,SF,IF.E);


	NF = quadOp(FF,FF,[]);
%	SF = (NF(1)-NF(2)+NF(3)-NF(4))/2; %Further jobs
%	LF = (NF(1)+NF(2)-NF(3)-NF(4))/2; %

	path1 = ['/work/dinnercha/NRG_data/','j',getenv('SLURM_ARRAY_JOB_ID'),'t',getenv('SLURM_ARRAY_TASK_ID'),'test']
	path2 = go('data')
    nrgdata = path1
	%nrgdata = []%memorymode?
	HU = QSpace;
	Hepsd = QSpace; %Make empty QSpace

	for i = 1:numel(FF)
		HU = HU + (U/2)*(sum(NF(i)) * (sum(NF(i)) - EF)); %EF =IDentity

		Hepsd =  Hepsd + mu * NF(i); % Quadratic term i think.

	end

	A0 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]);

    HU = contract(A0,'!2*',{HU,'!1',A0});
    Hepsd = contract(A0,'!2*',{Hepsd,'!1',A0});

	

	
	for i = 1:iter

		disp("=======================")
		disp("Current step is => "+ i)
		disp("=======================")

		[ocont,SE,nrgdataz,Aconts,ImSEle] = impSE2(Hepsd,HU,A0,FF,ZF,RhoV2in_arr(:,:,end),T,Lambda,nz,Nkeep,nrgdata,'numCorr',numCorr,'doZLD',{'minTE',0.0,'NaNtol',1e-40,'Nfit',inf});
		[nd,dd] = getEpVal(nrgdataz,sum(NF),sum(NF)*(sum(NF)-EF)/2)

		[RhoV2,Delta,A,G] = Bethe(ocont,mu,'SE',SE);
	
		G_arr = cat(3,G_arr,G);
		A_arr = cat(3,A_arr,A);
		SE_arr = cat(3,SE_arr,SE);
		Delata_arr = cat(3,Delta_arr,Delta);
		RhoV2_arr = cat(3,RhoV2_arr,RhoV2);
		nd_arr(end+1) = nd
		dd_arr(end+1) = dd
		
		[iscvg,RhoV2in] = updateHyb(RhoV2in_arr,RhoV2_arr,3,'Broyden',3,'-v');
		
		RhoV2in_arr = cat(3,RhoV2in_arr,RhoV2in);

		save(resname,'-v7.3')
		if iscvg

			disp("===================================================")
			disp("The DMFT calculation is converged at step =>  " + i)
			disp("================================================DONE")
			break
			
		end
		

	end
	%% energy flow

	[Etot,Qtot] = plotE(nrgdataz{end},'noshow');
	[ocont,SusConts] = impSusc(nrgdataz,[NF],T/5);%charge
	%[ocont,SusConts] = impSusc(nrgdataz,[NF;SF;LF],T/5);%charge, Spin, orbit
	
	

	save(resname,'-v7.3')
	disp(resname)
	disp('SAVE DONE!!')

	%%Susceptibility 
	%% spin, charge, 
	%% 
	%% 다른 quantity
end
