function benchmark_par (varargin)
    parfn = 'TEST'; % Parameterfile name 
    while numel(varargin) > 0
        if ischar(varargin{1})
            parfn = varargin{1};
            varargin(1) = [];
        else
            disp(varargin{1});
            error('ERR: Unknown input/option.');
        end
    end

    PE = 4; % # of cores to be occupied in clusters
    h_vmem = 64; % Memory (in GB) to be occupied in clusters
    syms = cell(1,0); % non-Abelian symmetry types to be exploited


    U = (0:1e-3:1e-2); % interaction strength
    % epsd = -U / 2; % on-site level
    % B = 0; % Zeeman field at the impurity, along the z direction

    T = 1e-8; % Temperature

    % Hybridization function parametrized by the frequency grid 'ozin' and the
    % function value 'RhoV2in' evaluated at 'ozin'. Here consider a simple
    % box-shaped case.
    D = 1; % half-bandwidth
    Gamma = 1e-4; % hybridization strength
    ozin = [-D; D];
    RhoV2in = (Gamma / pi) * [1; 1]; % values outside of the 'ozin' grid are assumed to be zero

    % NRG parameter
    Lambda = 2;

    partot = struct;
    cnt = 0;
    for it1 = (1:numel(U))
            cnt = cnt + 1;
            partot(cnt).U = U(it1)
            partot(cnt).T = T;
            partot(cnt).Lambda = Lambda;
            partot(cnt).RhoV2in = RhoV2in;
    end
    partot = partot(:);

    dispstruct(partot);
    disp(partot(1))

    parfn = [go('mu/Para/'),parfn,'.mat'];
    save(parfn,'partot','h_vmem','PE','syms','-v7.3');
    disp(['Saved to : ',parfn]);
end
