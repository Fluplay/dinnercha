function test_par (varargin)
    parfn = '2orb_Bethe_test_nz2'; % Parameterfile name ?
    while numel(varargin) > 0
        if ischar(varargin{1})
            parfn = varargin{1};
            varargin(1) = [];
        else
            disp(varargin{1});
            error('ERR: Unknown input/option.');
        end
    end

    PE = 8; % # of cores to be occupied in clusters
    h_vmem = 48; % Memory (in GB) to be occupied in clusters
    syms = cell(1,0); % non-Abelian symmetry types to be exploited
	D = 2 %Half Band Width
    U = (0.0:1:10); %D is setted to be 2 ; 0D to 4D
    mu = [-U/2.0 -U/2.0];
    T = 1e-8;
    Lambda = 4;
    Nkeep = 1600;
    nz = 2;
	
	ocont = getAcont(0,0,0,0)
	[ocont,RhoV2init] = get_Bethe(ocont,D)
    %RhoV2init = [];
	RhoV2init = RhoV2init'
	RhoV2init_arr = cat(2,RhoV2init,RhoV2init)
    partot = struct;
    cnt = 0;
    for it1 = (1:numel(U))
            cnt = cnt + 1;
            partot(cnt).U = U(it1);
            partot(cnt).mu = [-U(it1)/2.0 -U(it1)/2.0];
            partot(cnt).T = T;
            partot(cnt).Lambda = Lambda;
            partot(cnt).Nkeep = Nkeep;
            partot(cnt).nz = nz;
            partot(cnt).RhoV2init = RhoV2init_arr;
			partot(cnt).ocont = ocont
    end
    partot = partot(:);

    dispstruct(partot);
    disp(partot(1))

    parfn = [go('mu/Para/'),parfn,'.mat'];
    save(parfn,'partot','h_vmem','PE','syms','-v7.3');
    disp(['Saved to : ',parfn]);
end
