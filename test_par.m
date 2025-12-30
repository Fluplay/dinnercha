function test_par (varargin)
    parfn = 'test';
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

    U = (0.1:0.5:1);
    mu = U/2.0;
    T = exp((-10:5:0));
    Lambda = 4;
    Nkeep = 1000;
    nz = 2;
    RhoV2init = [];

    partot = struct;
    cnt = 0;
    for it1 = (1:numel(U))
        for it2 = (1:numel(T))
            cnt = cnt + 1;
            partot(cnt).U = U(it1);
            partot(cnt).mu = mu(it1);
            partot(cnt).T = T(it2);
            partot(cnt).Lambda = Lambda;
            partot(cnt).Nkeep = Nkeep;
            partot(cnt).nz = nz;
            partot(cnt).RhoV2init = RhoV2init;
        end
    end
    partot = partot(:);

    dispstruct(partot);
    disp(partot(1))

    parfn = [go('mu/Para/'),parfn,'.mat'];
    save(parfn,'partot','h_vmem','PE','syms','-v7.3');
    disp(['Saved to : ',parfn]);
end
