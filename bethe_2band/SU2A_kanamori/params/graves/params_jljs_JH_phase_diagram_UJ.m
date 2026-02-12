function Nk4k_JH_phase_UJ (varargin)
	% have init, alpha 0 ver
    parfn = 'Nk4k_JH_phase_UJ'; % Parameterfile name ?
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
	D = 1 %Half Band Width
    U = (2.2:0.2:4); %D is setted to be 2 ; 0D to 4D
    %mu = [-U/2.0 -U/2.0];
    T = 1e-8;
    Lambda = 6;%
    Nkeep = 4000;%
    nz = 2;%

	alphaL = 0;

	nd = 2;


	Jl = 1;
	Js = 1;

	Ju = (0:0.025:0.25);
	D = 1;

	ocont = getAcont(0,0,0,0);
%	load('RV2s/init_zeeman.mat');
	[ocont,RhoV2init] = get_Bethe(ocont,D);
	RhoV2init = RhoV2init';
	RhoV2init_arr = cat(2,RhoV2init,RhoV2init);
	disp('comp')



    partot = struct;
    cnt = 0;
    for it1 = (1:numel(U))
		for it2 = (1:numel(Ju))
			cnt = cnt + 1;
			partot(cnt).U = U(it1);
			partot(cnt).Jl = U(it1)*Ju(it2);
			partot(cnt).Js = U(it1)*Ju(it2);
			partot(cnt).alphaL = alphaL;
			partot(cnt).nd = nd;

			partot(cnt).mu = 1e-40;
			partot(cnt).T = T;
			partot(cnt).Lambda = Lambda;
			partot(cnt).Nkeep = Nkeep;
			partot(cnt).nz = nz;
			partot(cnt).RhoV2init = RhoV2init_arr;
			partot(cnt).ocont = ocont;
		end
    end

    partot = partot(:);

    dispstruct(partot);
    disp(partot(1))

    parfn = [go('mu/Para/'),parfn,'.mat'];
    save(parfn,'partot','h_vmem','PE','syms','-v7.3');
    disp(['Saved to : ',parfn]);
end
