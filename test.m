function test (parfn, varargin)
try 
    setenv('RC_STORE',go('rcs')); % use local directory only
    partot = job_func_preamble(parfn,varargin{:});
    for it = (1:numel(partot))
        [ U , mu , T , Lambda , Nkeep , nz,  RhoV2init ] = loadvar(partot(it), ...
        {'U','mu','T','Lambda','Nkeep','nz','RhoV2init'}, ...
        { [], [] , [], []     , []    , [] , []        });

        % % result file path
        resname = go(['glo/DMFT/test_',par2str('U',U,'mu',mu,'T',T,'Lambda',Lambda,'Nkeep',Nkeep,'nz',nz)]);
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
         
        if ~isempty(RhoV2init)
            RhoV2ins(:,:,:,:,1) = interp1(RhoV2init{1},RhoV2init{2},ocont,'linear','extrap');
        elseif exist(resname,'file')
            disp2('Previous result from the same job is detected. Continue from its last result.');
            S = load3(resname,'RhoV2ins','itd','-v');
            RhoV2ins(:,:,1) = S.RhoV2ins(:,:,S.itd);
        end

        %%%
        disp(U);


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
