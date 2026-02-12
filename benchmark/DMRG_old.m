function DMRG_old(parfn, varargin)
try 
    setenv('RC_STORE',go('rcs')); % use local directory only
    partot = job_func_preamble(parfn,varargin{:});
    for it = (1:numel(partot))
        [Nkeep ] = loadvar(partot(it), ...
        {'Nkeep'}, ...
        { []});

        % % result file path
        resname = go(['glo/DMFT/Bethe_single/test','DMRG']);

		resname_DMFT = resname; %Resname for DMFT JH

        if ~isempty(getenv('SLURM_JOB_ID'))
            resname = [resname,'_j',getenv('SLURM_ARRAY_JOB_ID'),'_',getenv('SLURM_ARRAY_TASK_ID'),'.mat'];
        else
            resname = [resname,'_',datestr(now,'yy.mm.dd_HH-MM-SS'),'.mat'];
        end

        %%%
		tic
		Nsite = 40;
		J = 1;

		MPO = XYMPO(J, Nsite);

		[Minit, Egs, Eeig, Ulist, Hnowlist, Hmatout] = MPSinitialize(MPO, Nkeep, 'E', 'Enum', 1);


		Nsweep = 50;
		[MDMRG,E0,Eiter,Sv] = DMRG_GS_1site_QSpace (Minit,MPO,Nkeep,Nsweep);
		toc

		%%%
        save(resname,'-v7.3');
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
