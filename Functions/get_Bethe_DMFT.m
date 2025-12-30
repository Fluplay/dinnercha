function get_Bethe_DMFT(U,mu,N,iter,numCorr,D,resname,RhoV2init,T,Lambda,Nkeep,nz)

	%% todo
	%% Numcorr => should be optional
	%% the input for filling factor should beadjusted
	%%% => ex ; 'n', 0.5 => mu = -U/2
	% I think 18-32 lines can be out from the for loop; let's try
	% THIS IS PROTOTYPE 

	%%todo end


	[FF,ZF,SF,IF] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1)
		
	[FF,ZF,SF,EF] = setItag('s00','op',FF(:),ZF,SF(:),IF.E);
	
	NF = quadOp(FF,FF,[]);

	HU = (U/2) * sum(NF) * (sum(NF)-EF);

	Hmu = sum(NF) * mu;

	vac = setItag('L00',getvac(EF)) ;%Vaccum Identity

	A0 = getIdentity(vac,2,EF,2,'K00*',[1,3,2]);

	HU = contract(A0,'!2*',{HU,'!1',A0});
	Hmu = contract(A0,'!2*',{Hmu,'!1',A0});
	

    nrgdata = go('data/mytest'); %Work folder
	for i = (1:iter)
		disp("=======================")
		disp("Current step is => "+ i)
		disp("=======================")
		resname_loop = [resname,'step_',num2str(i),'.mat']

		[ocont,SE,nrgdataz,Aconts,ImSEle] = impSE2(Hmu,HU,A0,FF,ZF,RhoV2init,T,Lambda,nz,Nkeep,nrgdata,'numCorr',numCorr,'Nchain',76,'doZLD',{'NaNtol',1e-40,'Nfit',inf});

		[RhoV2,Delta,A,G] = Bethe(ocont,mu,'SE',SE,'D',D);
		
		RhoV2init = RhoV2;

		save(resname_loop,'-v7.3');
		

		if i == 1
			disp('first step does not show convergenceness')

		else

			Delta_diff = abs(Delta - Delta_old)
			meanDd = mean(Delta_diff)
			maxDd = max(Delta_diff)
			disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			disp('The differences btw old and new Delta is ')
			disp("|Delta(w) - Delta_old|")
			disp("Average => " + meanDd)
			disp("maximum => " + maxDd)
			disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			disp('The differences btw old and new Delta is ')

		end

		Delta_old = Delta;
	end
	
end
