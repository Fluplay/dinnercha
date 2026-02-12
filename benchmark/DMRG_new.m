function DMRG_new(parfn, varargin)
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
		monitor_memory('before')
		[MDMRG,E0,Eiter,Sv] = DMRG_GS_1site_QSpace (Minit,MPO,Nkeep,Nsweep);
		toc
		monitor_memory('after')

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
function monitor_memory(label)
    % 사용법: monitor_memory('Step 1 시작');

    % 1. 내 프로세스 ID(PID) 확인
    pid = feature('getpid');

    % 2. 리눅스 시스템 정보(/proc) 읽기
    % VmRSS: 현재 실제 물리 메모리 점유량
    % VmPeak: 프로그램 시작 후 찍은 최대 메모리량
    cmd = sprintf('grep -E "VmRSS|VmPeak" /proc/%d/status', pid);
    [status, result] = system(cmd);

    if status == 0
        % 결과 파싱 (KB -> GB 변환)
        lines = splitlines(strtrim(result));

        peak_str = extractAfter(lines{contains(lines, 'VmPeak')}, ':');
        rss_str  = extractAfter(lines{contains(lines, 'VmRSS')}, ':');

        peak_gb = sscanf(peak_str, '%f') / (1024^2);
        rss_gb  = sscanf(rss_str, '%f')  / (1024^2);

        fprintf('[PID: %d | %s] 현재(RSS): %.2f GB / 최대(Peak): %.2f GB\n', ...
                pid, label, rss_gb, peak_gb);
    else
        fprintf('메모리 정보를 읽을 수 없습니다 (Linux 환경 아님?).\n');
    end
end
