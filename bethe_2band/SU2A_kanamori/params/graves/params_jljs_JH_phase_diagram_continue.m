function Nk4k_JH_phase (varargin)
	% have init, alpha 0 ver
    parfn = 'Nk4k_JH_phase_patch3'; % Parameterfile name ?
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
    h_vmem = 64; % Memory (in GB) to be occupied in clusters
    syms = cell(1,0); % non-Abelian symmetry types to be exploited
	D = 1 %Half Band Width
    %U = (0.1:0.2:2); %D is setted to be 2 ; 0D to 4D
    U = 1.1; %
    %mu = [-U/2.0 -U/2.0];
    T = 1e-8;
    Lambda = 6;%
    Nkeep = 4000;%
    nz = 2;%

	alphaL = 1;

	nd = 2;


	Jl = 1;
	Js = 1;

	Ju = [0.145,0.155,0.16,0.17,0.175];
	D = 1;

	ocont = getAcont(0,0,0,0);
	load('RV2s/init_zeeman.mat');
	[ocont,RhoV2init] = get_Bethe(ocont,D);
	RhoV2init = RhoV2init';
	RhoV2init_arr = cat(2,RhoV2init,RhoV2init);
	disp('comp')

	
	paths1 = 'Bethe_jljs_U=1.1_Jl=0.1595_Js=0.1595_T=1e-08_Lambda=6_Nkeep=4000_nz=2_alphaL=1_nd=2_D=_j155231_10.mat'
	paths2 = 'Bethe_jljs_U=1.1_Jl=0.1705_Js=0.1705_T=1e-08_Lambda=6_Nkeep=4000_nz=2_alphaL=1_nd=2_D=_j154829_2.mat'
	paths3 = 'Bethe_jljs_U=1.1_Jl=0.176_Js=0.176_T=1e-08_Lambda=6_Nkeep=4000_nz=2_alphaL=1_nd=2_D=_j154829_3.mat'
	paths4 = 'Bethe_jljs_U=1.1_Jl=0.187_Js=0.187_T=1e-08_Lambda=6_Nkeep=4000_nz=2_alphaL=1_nd=2_D=_j154829_5.mat'
	paths5 = 'Bethe_jljs_U=1.1_Jl=0.1925_Js=0.1925_T=1e-08_Lambda=6_Nkeep=4000_nz=2_alphaL=1_nd=2_D=_j154829_6.mat'

	prev1 = load(paths1)
	prev2 = load(paths2)
	prev3 = load(paths3)
	prev4 = load(paths4)
	prev5 = load(paths5)

    partot = struct;
    cnt = 0;
    for it1 = (1:5)
			cnt = cnt + 1;
			partot(cnt).U = U;
			partot(cnt).Jl = U*Ju(it1);
			partot(cnt).Js = U*Ju(it1);
			partot(cnt).alphaL = alphaL;
			partot(cnt).nd = nd;

			partot(cnt).mu = 1e-40;
			partot(cnt).T = T;
			partot(cnt).Lambda = Lambda;
			partot(cnt).Nkeep = Nkeep;
			partot(cnt).nz = nz;
			partot(cnt).RhoV2init = RhoV2init_arr;
			partot(cnt).ocont = ocont;

			if it1 == 1
				partot(cnt).prev = prev1
			end 

			if it1 == 2
				partot(cnt).prev = prev2
			end 

			if it1 == 3
				partot(cnt).prev = prev3
			end 

			if it1 == 4
				partot(cnt).prev = prev4
			end

			if it1 == 5
				partot(cnt).prev = prev5
			end
				
    end

    partot = partot(:);

    dispstruct(partot);
    disp(partot(1))

    parfn = [go('mu/Para/'),parfn,'.mat'];
    save(parfn,'partot','h_vmem','PE','syms','-v7.3');
    disp(['Saved to : ',parfn]);
end
