function mv(varargin)
    % MV - Changes the current directory to the path of the specified file
    %
    % Syntax:
    %   mv filePath
    %
    % Input:
    %   filePath - Full path to a file (without parentheses)
    %
    % Description:
    %   This function extracts the directory path from the given file path
    %   and changes the current directory to that location.
    %   It can be called without parentheses in MATLAB command window.
    %
    % Example:
    %   mv /home/subini0213/MATLABuser/Workspace/orbitalHeisenberg/compiled/JKI_VUMPS.m
    %   % This will change the current directory to:
    %   % /home/subini0213/MATLABuser/Workspace/orbitalHeisenberg/compiled/
    
        % Check if input is provided
        if nargin < 1
            error('File path must be provided');
        end
        
        % Combine all input arguments to handle spaces in path
        if iscell(varargin)
            filePath = strjoin(varargin, ' ');
        else
            filePath = varargin;
        end
        
        % Get the directory path by removing the file name
        [dirPath, ~, ~] = fileparts(filePath);
        
        % Check if the directory exists
        if ~exist(dirPath, 'dir')
            error('Directory "%s" does not exist', dirPath);
        end
        
        % Change directory
        cd(dirPath);
        
        % Display confirmation message
        fprintf('Current directory changed to: %s\n', dirPath);
    end