function get_Bethe_DMFT_jljs_AA(U, ...
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
	path1 = ['/work/dinnercha/NRG_data/','j',getenv('SLURM_ARRAY_JOB_ID'),'t',getenv('SLURM_ARRAY_TASK_ID'),'test'];
	nrgdata = path1
    RhoV2_arr = [];
    RhoV2in_arr = RhoV2init;
    G_arr = [];
    A_arr = [];
    SE_arr = [];
    Delta_arr = [];
    nd_arr = [];
    dd_arr = [];

	[FF,ZF,~,IF] = getLocalSpace('FermionS','Acharge(:),Aspin(:)','NC',2);
	[FF,ZF,EF] = setItag('s00','op',FF,ZF,IF.E);
	NF = quadOp(FF,FF,[]);
	Ntot = sum(NF(:)); % #tot op

    Sx = (quadOp(FF(1,1),FF(1,2),'*')+quadOp(FF(1,2),FF(1,1),'*')+quadOp(FF(2,1),FF(2,2),'*')+quadOp(FF(2,2),FF(2,1),'*'))/2;

    Sy = (quadOp(FF(1,2),FF(1,1),'*')-quadOp(FF(1,1),FF(1,2),'*')+quadOp(FF(2,2),FF(2,1),'*')-quadOp(FF(2,1),FF(2,2),'*'))*(1i/2);

    SF = (NF(1,1)-NF(1,2)+NF(2,1)-NF(2,2))/2; % == SF(3); Sz?
    SF2 = quadOp(Sx,Sx,[])+quadOp(Sy,Sy,[])+quadOp(SF,SF,[]);

    LF = (NF(1,1)+NF(1,2)-NF(2,1)-NF(2,2))/2;
	Lx = (quadOp(FF(1,1),FF(2,1),'*') + quadOp(FF(1,2),FF(2,2),'*') + quadOp(FF(2,1),FF(1,1),'*') + quadOp(FF(2,2),FF(1,2),'*'))/2
	Ly = ( (quadOp(FF(2,1),FF(1,1),'*') + quadOp(FF(2,2),FF(1,2),'*')) - ...
       (quadOp(FF(1,1),FF(2,1),'*') + quadOp(FF(1,2),FF(2,2),'*')) ) * (1i/2)	

	L2 = quadOp(Lx,Lx,[]) + quadOp(Ly,Ly,[]) + quadOp(LF,LF,[])

	
    HU = (U-3/2*Js)/2*(Ntot-nd*EF)*(Ntot-nd*EF) + Js*(SF2) - Jl*(L2 + alphaL*LF*LF);

	Hepsd = -mu*Ntot


	A0 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]);

	HU = contract(A0,'!2*',{HU,'!1',A0});
	Hepsd = contract(A0,'!2*',{Hepsd,'!1',A0});
	
	FF_input = QSpace()
	FF_input(1) = FF(1,1)
	FF_input(2) = FF(1,2)
	FF_input(3) = FF(2,1)
	FF_input(4) = FF(2,2)

	for i = 1:iter

        disp("=======================")
        disp("Current step is => "+ i)
        disp("=======================")

        [ocont,SE,nrgdataz,Aconts,ImSEle] = impSE2(Hepsd,HU,A0,FF_input,ZF,RhoV2in_arr(:,:,end),T,Lambda,nz,Nkeep,nrgdata,'numCorr',numCorr,'doZLD',{'minTE',0.0,'NaNtol',1e-40,'Nfit',inf});
        [nd,dd] = getEpVal(nrgdataz,NF,sum(NF(:))*(sum(NF(:))-EF)/2);
		disp("nd idx : (orbit,spin)")
		disp("Expected nd =>" + nd)
		disp("Expected dd =>" + dd)
        [RhoV2,Delta,A,G] = Bethe(ocont,mu,'SE',SE);

        G_arr = cat(3,G_arr,G);
        A_arr = cat(3,A_arr,A);
        SE_arr = cat(3,SE_arr,SE);
        Delata_arr = cat(3,Delta_arr,Delta);
        RhoV2_arr = cat(3,RhoV2_arr,RhoV2);
		nd_arr = cat(3,nd_arr,nd)
        dd_arr(end+1) = dd

        [iscvg,RhoV2in] = updateHyb(RhoV2in_arr,RhoV2_arr,3,'Broyden',5,'-v');

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
    [ocont,SusConts_NF11] = impSusc(nrgdataz,[NF(1,1)],T/5);%charge
    [ocont,SusConts_NF12] = impSusc(nrgdataz,[NF(1,2)],T/5);%charge
    [ocont,SusConts_NF21] = impSusc(nrgdataz,[NF(2,1)],T/5);%charge
    [ocont,SusConts_NF22] = impSusc(nrgdataz,[NF(2,2)],T/5);%charge
	




    save(resname,'-v7.3')
    disp(resname)
    disp('SAVE DONE!!')
end	
