function GDZ_TESTU3 (varargin)
	% have init, alpha 0 ver
    parfn = 'GDZ_TESTU3'; % Parameterfile name ?
    while numel(varargin) > 0
        if ischar(varargin{1})
            parfn = varargin{1};
            varargin(1) = [];
        else
            disp(varargin{1});
            error('ERR: Unknown input/option.');
        end
    end

    PE = 6; % # of cores to be occupied in clusters
    h_vmem = 64; % Memory (in GB) to be occupied in clusters
    syms = cell(1,0); % non-Abelian symmetry types to be exploited
	D = 1 %Half Band Width
    U = 3; %D is setted to be 2 ; 0D to 4D
    %mu = [-U/2.0 -U/2.0];
    T = 1e-8;
    Lambda = 4;%
    Nkeep = 4000;%
    nz = 2;%

	alphaL = 0;



	%Jd = 1;
	%Js = 1;

	Ju = (-0.2:0.05:0.2); % case for Jd = Js
	D = 1;

	ocont = getAcont(0,0,0,0);
	[ocont,RhoV2init] = get_Bethe(ocont,D);
	RhoV2init = RhoV2init';
	RhoV2init_arr = cat(2,RhoV2init,RhoV2init);

    partot = struct;
    cnt = 0;
    for it1 = (1:numel(Ju))
		Jd = U*Ju(it1);
		Js = U*Ju(it1);
		cnt = cnt + 1;
		partot(cnt).U = U;
		partot(cnt).Jd = Jd;
		partot(cnt).Js = Js;
		partot(cnt).alphaL = alphaL;


		partot(cnt).mu = -3/2*U+1/4*Js+1/2*Jd;
		partot(cnt).T = T;
		partot(cnt).Lambda = Lambda;
		partot(cnt).Nkeep = Nkeep;
		partot(cnt).nz = nz;
		partot(cnt).RhoV2init = RhoV2init_arr;
		partot(cnt).ocont = ocont;
    end

    partot = partot(:);

    dispstruct(partot);
    disp(partot(1))

    parfn = [go('mu/Para/'),parfn,'.mat'];
    save(parfn,'partot','h_vmem','PE','syms','-v7.3');
    disp(['Saved to : ',parfn]);
end
