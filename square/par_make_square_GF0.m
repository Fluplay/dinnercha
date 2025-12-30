function GF0_par (varargin)
    parfn = 'make_Sqaure_GF0'; % Parameterfile name ?
    while numel(varargin) > 0
        if ischar(varargin{1})
            parfn = varargin{1};
            varargin(1) = [];
        else
            disp(varargin{1});
            error('ERR: Unknown input/option.');
        end
    end

    PE = 32; % # of cores to be occupied in clusters
    h_vmem = 256; % Memory (in GB) to be occupied in clusters
    syms = cell(1,0); % non-Abelian symmetry types to be exploited
	D = 2 %Half Band Width
    U = 0.0; %D is setted to be 2 ; 0D to 4D
    mu = -U/2.0;
    T = 1e-8;
    Lambda = 4;
    Nkeep = 2000;
    nz = 4;
	
	ocont = getAcont(0,0,0,0)
	[ocont,RhoV2init] = get_Bethe(ocont,D)
    %RhoV2init = [];

    partot = struct;
   	partot.U = U
    partot.mu = mu;
    partot.T = T;
	partot.Lambda = Lambda;
    partot.Nkeep = Nkeep;
    partot.nz = nz;
    partot.RhoV2init = RhoV2init;
	partot.ocont = ocont
    partot = partot(:);

    dispstruct(partot);
    disp(partot(1))

    parfn = [go('mu/Para/'),parfn,'.mat'];
    save(parfn,'partot','h_vmem','PE','syms','-v7.3');
    disp(['Saved to : ',parfn]);
end
