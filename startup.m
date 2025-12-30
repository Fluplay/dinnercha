%% Startup file that sets up environment variables and add paths
% Revise the blocks enclosed between "TODO (Start)" and "TODO (End)", to be
% suitable for an individual user's situation.

%% Set environment variables
% Below there are several environment variables (starting with $) that
% specify the paths to directories used in MuNRG. To understand them and
% set your own, you should be familiar with the file systems of different
% cluster systems you use.
% Once you set the paths correctly, you can use 'go' function for changing
% to directories specified by those variables; see below and the
% documentation of 'go' for details.

% $LOC_CLUSTER: A shorthand name for a system. When you work on multiple
%       systems (e.g., your laptop, local cluster in your department, and a
%       supercomputer), you can set distinct values of $LOC_CLUSTER based
%       on their OS and hostnames. $LOC_CLUSTER can be used to
%       differentiate paths to directories used in MuNRG (see below).
%       (NOTE: Some cluster submission related functions recognize some
%       preset values of $LOC_CLUSTER ('ASC', 'LRZ', or 'SNU'). We plan to
%       clean them up so that other users can better adapt them to their
%       systems.)
%%%%% TODO (Start) %%%%%
% obtain the static hostname (just querying 'hostname' can give varying
% results, e.g., if one uses a WiFi hotspot.)
if strcmp(computer,'MACI64')
    hname = strtrim(evalc('system(''scutil --get ComputerName'');'));
else
    hname = strtrim(evalc('system(''hostname -f'');'));
end

if strcmp(computer,'MACI64') && strcmp(hname,'Seung-Sup’s MacBook Pro')
    setenv('LOC_CLUSTER','HOME'); % S.Lee's MacBook
    % Note: If you use the Apple silicon version of MATLAB, use
    % 'strcmp(computer,'MACA64')'
elseif strcmp(computer,'GLNXA64') && (strcmp(hname,'headnode') || ...
        ~isempty(regexp(hname,'^[a-z]\d\d$','once')))
    % if it's a Linux machine with hostname 'headnode' or (regexp) '[a-z]\d\d' (e.g. b21)
    setenv('LOC_CLUSTER','SNU'); % SNU CMP theory cluster
elseif strcmp(computer,'PCWIN64') && strcmp(hname,'SEUNGSUPLEEFA4D')
    setenv('LOC_CLUSTER','HOME_WIN'); % S.Lee's Windows running on a MacBook
else
    error('ERR: Unknown setting?')
end
%%%%% TODO (End) %%%%%

% $MUNRG_DIR: the path for the MuNRG directory; to be referred by go('mu').
setenv('MUNRG_DIR',fileparts(mfilename('fullpath'))); % directory that contains this startup.m

% $QSPACE_DIR: the path for the QSpace directory; to be referred by
%       go('qs'). The example below assumes the situation that the QSpace 
%       and MuNRG directories are the subdirectories of the same directory
%       with the same depth level.
%%%%% TODO (Start) %%%%%
setenv('QSPACE_DIR',[fileparts(getenv('MUNRG_DIR')),filesep,'QSpace_v4']); % QSpace path
% NOTE: Don't put '/' at the end of the path.
%%%%% TODO (End) %%%%%

% $GLOBAL_DIR: The directory to be accessed by all computational nodes for
%       file I/O. It contains input parameter files and the central
%       Clebsch-Gordan coefficient tensor (CGT) database. Referred by
%       go('glo') or go('pts').
% $SCRATCH_DIR: The path for a directory containing temporary files (e.g.
%       raw NRG data, job-dependently created CGT data) used in a single
%       session. In case of cluster jobs, $SCRATCH_DIR is set differently
%       for all different (sub-)tasks to avoid interference. There are
%       cluster systems where nodes don't have local scratch drives and a
%       network drive for global scratch (to be specified by $GLOBAL_DIR)
%       exists instead. For such systems, $SCRATCH_DIR can be a
%       sub-directory of $GLOBAL_DIR. Referred by go('tmp') or go('data').
% $STORAGE_DIR: The path to store calculation results for a longer period.
%       Referred by go('sto') or go('pt').
% NOTE: Don't put '/' at the end of those paths.
%%%%% TODO (Start) %%%%%
switch getenv('LOC_CLUSTER')
    case 'HOME'
        setenv('GLOBAL_DIR',[getenv('HOME'),filesep,'data']); % use a directory under the home directory
        setenv('SCRATCH_DIR',[getenv('GLOBAL_DIR'),filesep,'scratch']);
        setenv('STORAGE_DIR',getenv('GLOBAL_DIR'));
    case 'SNU'
        setenv('GLOBAL_DIR',['/project/',getenv('USER')]);
        setenv('SCRATCH_DIR',['/tmp/',getenv('USER')]);
        setenv('STORAGE_DIR',['/data/',getenv('USER')]);
    case 'HOME_WIN'
        setenv('GLOBAL_DIR',['Y:',filesep,'data']); % Mac's home is linked as the Y drive
        setenv('SCRATCH_DIR',[getenv('GLOBAL_DIR'),filesep,'scratch']);
        setenv('STORAGE_DIR',getenv('GLOBAL_DIR'));
    otherwise
        error('ERR: Unknown setting?')
end
%%%%% TODO (End) %%%%%

% If running on a cluster via SLURM, use separate temporary directories and
% log files for different array jobs. Those directories are named as e.g.,
% j1234t56.
if ~isempty(getenv('SLURM_JOB_ID')) || ~isempty(getenv('SLURM_ARRAY_JOB_ID')) % detects whether it's on SLURM
    if ~isempty(getenv('SLURM_ARRAY_TASK_ID'))
        setenv('SCRATCH_DIR',[getenv('SCRATCH_DIR'),filesep,'j',getenv('SLURM_ARRAY_JOB_ID'), ...
            't',getenv('SLURM_ARRAY_TASK_ID')]);
        setenv('LOG_FILE',[getenv('GLOBAL_DIR'),filesep,getenv('SLURM_JOB_NAME'),'_j',getenv('SLURM_ARRAY_JOB_ID'), ...
            't',getenv('SLURM_ARRAY_TASK_ID'),'.log']);
    else % for a system that does not allow array tasks, one can use 'mpiexec' to spawn them
        setenv('SCRATCH_DIR',[getenv('SCRATCH_DIR'),filesep,'j',getenv('SLURM_JOB_ID'), ...
            't',sprintf('%i',str2double(getenv('PMI_RANK'))+1)]);
        setenv('LOG_FILE',[getenv('GLOBAL_DIR'),filesep,getenv('SLURM_JOB_NAME'),'_j',getenv('SLURM_JOB_ID'), ...
            't',sprintf('%i',str2double(getenv('PMI_RANK'))+1),'.log']);
    end
end

% $RC_STORE: the path for the CGT data to be used by the QSpace library.
%       For clusters, this path contains local CGT data. Referred by
%       go('rcs').
setenv('RC_STORE',[getenv('SCRATCH_DIR'),filesep,'RCStore']);

% create directories
strs = {'GLOBAL_DIR','SCRATCH_DIR','STORAGE_DIR','RC_STORE'};
for its = (1:numel(strs))
    if ~exist(getenv(strs{its}),'dir')
        fprintf(['mkdir: ',getenv(strs{its}),'\n']);
        mkdir(getenv(strs{its}));
    end
end

% Set the path for the central CGT database. QSpace generates CGT data on
% the fly, so CGT data generated by parallel jobs may have conflicts. So
% MuNRG sets $RC_STORE as 'global_path:local_path', where the central CGT
% database is stored at 'global_path' and the job-dependent CGT data is at
% 'local_path'. For better performance, regularly update the central
% database with one of the local data; be careful about conflicts!
%%%%% TODO (Start) %%%%%
if ~strcmp(getenv('LOC_CLUSTER'),'HOME')
    setenv('RC_STORE',[getenv('GLOBAL_DIR'),'/RCSarxiv:',getenv('RC_STORE')]);
                      %         global_path           :    local_path
%%%%% TODO (End) %%%%%
end

% $MPNRG_BASE_PATH: The common base of the paths to be used in the
%       multipoint NRG calculations. By only changing this environment
%       variable, one can bring the data generated from one system (e.g.,
%       SuperMUC-NG) into another (e.g., ASC).
setenv('MPNRG_BASE_PATH',[getenv('GLOBAL_DIR'),filesep,'mpNRG']);

%% Set paths to be read by MATLAB
if ~(isdeployed || ismcc) % If the code runs within interactive MATLAB environment
    warning('off','MATLAB:dispatcher:pathWarning');
    restoredefaultpath;
    warning('on','MATLAB:dispatcher:pathWarning');
    
    % add the relevant sub-directories of QSpace to path
    strs = {'Class','bin','lib','util','setup','tensor'};
    for its = (1:numel(strs))
        addpath(genpath([getenv('QSPACE_DIR'),filesep,strs{its}]));
    end

    % add the MuNRG directory and its relevant sub-directories to path
    addpath(getenv('MUNRG_DIR'));
    strs = {'ClusterSubmissionScripts','Compiled','DMFT','Impurities',...
        'Multipoint','NRG','Para','Util'};
    for its = (1:numel(strs))
        addpath(genpath([getenv('MUNRG_DIR'),filesep,strs{its}]));
    end
    
    % You can also add another directories to path. For example, to put
    % your own codes using QSpace and MuNRG, you can create a
    % sub-directory, whose name is equal to your user name, under the MuNRG
    % directory:
    % addpath(genpath([getenv('MUNRG_DIR'),filesep,getenv('USER')]));
	%%%%% TODO (Start) %%%%%
    addpath(genpath([getenv('MUNRG_DIR'),filesep,'SLee']));
    %%%%% TODO (End) %%%%%
end

dispbox('-width',-88, ...
    '>> MuNRG/startup.m', ...
    [' Machine location: ',getenv('LOC_CLUSTER')], ...
    ['      QSpace path: ',getenv('QSPACE_DIR')], ...
    ['       MuNRG path: ',getenv('MUNRG_DIR')], ...
    ['  Global work dir: ',getenv('GLOBAL_DIR')], ...
    ['      Scratch dir: ',getenv('SCRATCH_DIR')], ...
    ['Long-term storage: ',getenv('STORAGE_DIR')], ...
    ['         CGT data: ',getenv('RC_STORE')]);

clear
