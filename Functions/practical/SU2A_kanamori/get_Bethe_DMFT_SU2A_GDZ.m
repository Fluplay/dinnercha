function get_Bethe_DMFT_SU2A_GDZ(U, ...
						Js, ...
						Jd, ...
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
						alphaL)
	path1 = ['/work/dinnercha/NRG_data/','j',getenv('SLURM_ARRAY_JOB_ID'),'t',getenv('SLURM_ARRAY_TASK_ID'),'test'];
	nrgdata = path1
    RhoV2_arr = [];
    RhoV2in_arr = RhoV2init;
    G_arr = [];
    A_arr = [];
    SE_arr = [];
    Delta_arr = [];
    dd_arr = [];
	nd_arr = [];
	Lz_arr = [];

	SF = 'Nothing'
	LF = 'Nothing'
	NF = 'Nothing'
	

	for i = 1:iter
		
		[Hepsd,HU,A0,FF,ZF,NF,LF,EF,SF] = get_ham(U,Jd,Js,alphaL,mu);

        disp("=======================")
        disp("Current step is => "+ i)
        disp("=======================")

        [ocont,SE,nrgdataz,Aconts,ImSEle] = impSE2(Hepsd,HU,A0,FF,ZF,RhoV2in_arr(:,:,end),T,Lambda,nz,Nkeep,nrgdata,'numCorr',numCorr,'doZLD',{'minTE',0.0,'NaNtol',1e-40,'Nfit',inf});
        [nn,dd] = getEpVal(nrgdataz,NF,sum(NF(:))*(sum(NF(:))-EF)/2);
		Lz = getEpVal(nrgdataz,LF)

		disp("nn idx : (orbit,spin)")
		disp("Expected nd =>" + nn)
		disp("Expected dd =>" + dd)

		%Symmetrize

		SEtmp = SE;
		SE(:,1) = (SEtmp(:,1) + SEtmp(:,2) )/2;
		SE(:,2) = SE(:,1);
		clear SEtmp


        [RhoV2,Delta,A,G] = Bethe(ocont,mu,'SE',SE);
			
        G_arr = cat(3,G_arr,G);
        A_arr = cat(3,A_arr,A);
        SE_arr = cat(3,SE_arr,SE);
        Delata_arr = cat(3,Delta_arr,Delta);
        RhoV2_arr = cat(3,RhoV2_arr,RhoV2);
		nd_arr = cat(3,nd_arr,nn)
		
        dd_arr(end+1) = dd
		Lz_arr(end+1) = Lz

        [iscvg,RhoV2in] = updateHyb(RhoV2in_arr,RhoV2_arr,3,'Broyden',5,'amix',0.8,'-v');

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
    [ocont,SusConts_NF] = impSusc(nrgdataz,[NF(1);NF(2)],T/5);%charge
  
    [ocont,SusConts_LF] = impSusc(nrgdataz,LF,T/5);%charge

    [ocont,SusConts_SF] = impSusc(nrgdataz,SF,T/5);%charge
	




    save(resname,'-v7.3')
    disp(resname)
    disp('SAVE DONE!!')
end	

function [Hepsd,HU,A0,FF,ZF,NF,LF,EF,SF] = get_ham (U,Jd,Js,alphaL,mu)

	[F, Z, S, I] = getLocalSpace('FermionS', 'Acharge,SU2spin');
	Fcell = {F, Z, S, I.E; F, Z, S, I.E};
	sym = {'A'};
	symnum = [2; 1];
	tagname = {'a1', 'a2'};

	for it1 = 1:size(Fcell, 2) %F,Z,S,I respectively

	  for it2 = 1:size(symnum, 1) %2,4 to F1, 1,3 to F2

		for it3 = 1:numel(sym) %A,SU2 respectively
		  Fcell{it2, it1} = addSymmetry(Fcell{it2, it1}, sym{it3}, 'pos', symnum(it2, it3));
		end

		Fcell{it2, it1} = setItag(tagname{it2}, 'op', Fcell{it2, it1});
	  end

	end

	A = getIdentity(Fcell{1, 4}, 2, Fcell{2, 4}, 2, 's00');
	FF = QSpace(1, 2); ZFs = QSpace(1, 2); SF = QSpace(1, 2); NF = QSpace(1, 2);
	FF(1) = contract(A, '!3*', {Fcell{1, 1}, '!1', {A, Fcell{2, 2}}}, [1 3 2]);
	FF(2) = contract(A, '!3*', {Fcell{2, 1}, '!1', A}, [1 3 2]);

	for it1 = 1:size(Fcell, 1)
	  ZFs(it1) = contract(A, '!3*', {Fcell{it1, 2}, '!1', A});
	  SF(it1) = contract(A, '!3*', {Fcell{it1, 3}, '!1', A}, [1 3 2]);
	  NF(it1) = quadOp(FF(it1), FF(it1), []);
	end
	EF = contract(A, '!3*', {Fcell{1, 4}, '!1', A});
	ZF = ZFs(1) * ZFs(2);


	NF = quadOp(FF,FF,[]);
	
	SF12 = quadOp(SF(1),SF(2),[])
	NF12 = quadOp(NF(1),NF(2),[])

	Ntot = sum(NF);

	LF = (NF(1)-NF(2))/2;

	HU = (U-Jd)/2*Ntot*(Ntot-1)+Js*SF12+(Jd-1/4*Js)*NF12

	Hepsd = mu*Ntot

    A0 = getIdentity(setItag('L00',getvac(EF)),2,EF,2,'K00*',[1 3 2]);

    HU = contract(A0,'!2*',{HU,'!1',A0});
    Hepsd = contract(A0,'!2*',{Hepsd,'!1',A0});
	
end
