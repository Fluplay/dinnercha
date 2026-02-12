function get_Bethe_DMFT_jljs_SU2A_Zeeman_nv(U, ...
						Js, ...
						Jl, ...
						mu, ...
						iter, ...
						numCorr, ...
						D, ...
						resname, ...
						RhoV2init, ...
						T, ...
						Lambda, ...
						Nkeep, ...
						nz, ...
						alphaL, ...
						nd)
	path1 = ['/project/dinnercha/NRG_data/','j',getenv('SLURM_ARRAY_JOB_ID'),'t',getenv('SLURM_ARRAY_TASK_ID'),'test'];
	nrgdata = path1
    RhoV2_arr = [];
    RhoV2in_arr = RhoV2init;
    G_arr = [];
    A_arr = [];
    SE_arr = [];
    Delta_arr = [];
    nd_arr = [];
    dd_arr = [];
	Lz_arr = [];
	Lsq_arr = [];

	
	hl = 1e-3

	for i = 1:iter
		
		if i > 30
			disp("hl is modified")
			hl = 0.0

		end

		[Hepsd,HU,A0,FF,ZF,NF,LF,L2,EF] = get_ham(U,Jl,Js,hl,nd,alphaL,mu)

        disp("=======================")
        disp("Current step is => "+ i)
        disp("Current hl is  => "+ hl)
        disp("=======================")

        [ocont,SE,nrgdataz,Aconts,ImSEle] = impSE2(Hepsd,HU,A0,FF,ZF,RhoV2in_arr(:,:,end),T,Lambda,nz,Nkeep,nrgdata,'numCorr',numCorr,'doZLD',{'minTE',0.0,'NaNtol',1e-40,'Nfit',inf});
        [nn,dd] = getEpVal(nrgdataz,NF,sum(NF(:))*(sum(NF(:))-EF)/2);
		Lz = getEpVal(nrgdataz,LF)
		Lsq = getEpVal(nrgdataz,L2)

		disp("nd idx : (orbit,spin)")
		disp("Expected nd =>" + nn)
		disp("Expected dd =>" + dd)
        [RhoV2,Delta,A,G] = Bethe(ocont,mu,'SE',SE);

        G_arr = cat(3,G_arr,G);
        A_arr = cat(3,A_arr,A);
        SE_arr = cat(3,SE_arr,SE);
        Delata_arr = cat(3,Delta_arr,Delta);
        RhoV2_arr = cat(3,RhoV2_arr,RhoV2);
		nd_arr = cat(3,nd_arr,nn)
		
        dd_arr(end+1) = dd
		Lz_arr(end+1) = Lz
		Lsq_arr(end+1) = Lsq

        [iscvg,RhoV2in] = updateHyb(RhoV2in_arr,RhoV2_arr,3,'Broyden',0,'amix',0.7,'-v');

        RhoV2in_arr = cat(3,RhoV2in_arr,RhoV2in);

        save(resname,'-v7.3')
        if iscvg

            disp("===================================================")
            disp("The DMFT calculation is converged at step =>  " + i)
            disp("================================================DONE")
            break

        end


    end

	disp("======Loop is Over. trying to calculate SuScept=====")
    [Etot,Qtot] = plotE(nrgdataz{end},'noshow');
    %[ocont,SusConts_SF] = impSusc(nrgdataz,[SF],T/5);%charge
    [ocont,SusConts_LFSF] = impSusc(nrgdataz,[LF;SF],T/5);%charge
	O = NF(:)

    [ocont,SusConts_LFSF] = impSusc(nrgdataz,[LF;SF;O],T/5);%charge
    [ocont,SusConts_NF11] = impSusc(nrgdataz,[NF(1)],T/5);%charge
    [ocont,SusConts_NF12] = impSusc(nrgdataz,[NF(2)],T/5);%charge
	




    save(resname,'-v7.3')
    disp(resname)
    disp('SAVE DONE!!')
end	

function [Hepsd,HU,A0,FF,ZF,NF,LF,L2,EF] = get_ham (U,Jl,Js,hl,nd,alphaL,mu)

	[FF,ZF,SF,IF] = getLocalSpace('FermionS','Acharge(:),SU2spin','NC',2);
	[FF,ZF,SF,EF] = setItag('s00','op',FF,ZF,SF,IF.E);
	NF = quadOp(FF,FF,[]);
	Ntot = sum(NF(:)); % #tot op

	NFsq = quadOp(NF,NF,[]);

	
	
	SF2 = quadOp(SF,SF,[]);

    LF = (NF(1) - NF(2) )/2;
	Lx = (quadOp(FF(1),FF(2),'*') + quadOp(FF(2),FF(1),'*'))/2;
	Ly =  (quadOp(FF(2),FF(1),'*') - quadOp(FF(1),FF(2),'*') ) * (1i/2)	;

	L2 = quadOp(Lx,Lx,[]) + quadOp(Ly,Ly,[]) + quadOp(LF,LF,[]);



	
    HU = (U-3/2*Js)/2*(Ntot-nd*EF)*(Ntot-nd*EF) + Js*(SF2) - Jl*(L2 + alphaL*LF*LF + hl*LF);
	
	Hepsd = -mu*Ntot;


	A0 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]);

	HU = contract(A0,'!2*',{HU,'!1',A0});
	Hepsd = contract(A0,'!2*',{Hepsd,'!1',A0});

end
