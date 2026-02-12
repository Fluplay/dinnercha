function SU2A_jljs_Zeeman_JH (parfn, varargin)
try
    setenv('RC_STORE',go('rcs')); % use local directory only
    partot = job_func_preamble(parfn,varargin{:});
    for it = (1:numel(partot))
        [ U ,Js, Jl , mu , T , Lambda , Nkeep , nz,  RhoV2init ocont, alphaL, nd, D] = loadvar(partot(it), ...
        {'U','Js','Jl','mu','T','Lambda','Nkeep','nz','RhoV2init','ocont','alphaL','nd','D'}, ...
        { [], [] , [], []     , [], [] , [], [],[],	 [], [],[] ,[]      });

        % % result file path
        %resname = go(['glo/DMFT/Bethe_single/Bethe_',par2str('U',U,'mu',mu,'T',T,'Lambda',Lambda,'Nkeep',Nkeep,'nz',nz)]);

        resname = go(['/local/MuNRG/mybin/bethe_2band/SU2A_kanamori/out/PD/Bethe_jljs_',par2str('U',U,'Jl',Jl,'Js',Js,'T',T,'Lambda',Lambda,'Nkeep',Nkeep,'nz',nz,'alphaL',alphaL,'nd',nd,'D',D)]);

		resname_DMFT = resname; %Resname for DMFT JH

        if ~isempty(getenv('SLURM_JOB_ID'))
            resname = [resname,'_j',getenv('SLURM_ARRAY_JOB_ID'),'_',getenv('SLURM_ARRAY_TASK_ID'),'.mat'];
        else
            resname = [resname,'_',datestr(now,'yy.mm.dd_HH-MM-SS'),'.mat'];
        end
        if ~exist(fileparts(resname),'dir')
            mkdir(fileparts(resname));
        end
        disp2(['Result will be saved to: ',resname]);

        nrgdata = cellfun(@(x) go(['data/NRG_itz=',sprintf('%i',x)]), num2cell(1:nz), 'UniformOutput', false);

        %if ~isempty(RhoV2init)
           % RhoV2ins(:,:,:,:,1) = interp1(RhoV2init{1},RhoV2init{2},ocont,'linear','extrap');
        %elseif exist(resname,'file')
        %    disp2('Previous result from the same job is detected. Continue from its last result.');
        %    S = load3(resname,'RhoV2ins','itd','-v');
        %    RhoV2ins(:,:,1) = S.RhoV2ins(:,:,S.itd);
        %end
		
		
		%%%
		disp("Comp")
		iter =30;
		numCorr = 3;
		D = 1;
		get_Bethe_DMFT_jljs_SU2A_Zeeman_init_JH(U,Js,Jl,mu,iter,numCorr,D,resname,RhoV2init,T,Lambda,Nkeep,nz,alphaL,nd)
		%%%




    end; clear it;
catch e
    disp2(getReport(e)); % Report error
    if ~(isdeployed || ismcc)
        keyboard
    end
    % Save current workspace
    errf = ['Error_',mfilename];
    if ~isempty(getenv('SLURM_ARRAY_JOB_ID'))
        errf = [errf,'_j',getenv('SLURM_ARRAY_JOB_ID')];
        if ~isempty(getenv('SLURM_ARRAY_TASK_ID'))
            errf = [errf,'t',getenv('SLURM_ARRAY_TASK_ID')];
        end
    end
    errf = [go('pts'),filesep,errf,'.mat'];
    save(errf,'-v7.3');
    dispbox('Error! Current workspace is saved in :',errf);
    % everything should be before rethrowing error
    rethrow(e);
end
end
