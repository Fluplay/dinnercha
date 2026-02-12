%% Test script for trying out NRG routines
clear;

num_threads_SL(8); % set this in the beginning; otherwise an error message occurs due to missing multithreading environment

%% Objective: Do NRG with SU(2)_c x SU(2)_s and identify RG fixed points

% % System paramters
% Impurity Hamiltonian parameters
N = max(ceil(-2 * log(T / 500) / log(Lambda)), 20);
nz = 3;
Nkeep = 300;
Nfit = round(-2 * log(1e-8) / log(Lambda));
Etrunc = []; % 9;
ETRUNC = []; % infs(1,20);

nrgdata = [];
% nrgdata = go('data/NRG/NRG'); % location to save the results in hard disk drive
[ff, gg, dff, dgg] = doZLD(ozin, RhoV2in, Lambda, N, nz, 'Nfit', Nfit);

[F, Z, S, I] = getLocalSpace('FermionS', 'SU2charge,SU2spin');

% Construct Hamiltonian
% SU2, SU2 symmetry
% Projectors into records with labels [1, 0], [0, 1]
Proj_10 = getsub(I.E, 1);
Proj_01 = getsub(I.E, 2);
assert(all(Proj_01.Q{1} == [0, 1])); % Check q-labels just in case

H0 = (-U / 2) * Proj_01; % On-site hamiltonian
H0.info.itags = {'K00', 'K00*'};

% Construct A matrix
% s00 is a vacuum space in this case
dummy = getvac(I.E);
A0 = permute(getIdentity(I.E, 1, dummy, 1), '231');
A0.info.itags = {'L00', 'K00*', 's00'};

% Call NRG routine
% Do not use gg in case of SU2charge
FC = {};
suffixes = transpose(reshape(sprintf('%02d', 0:N - 1), 2, []));
chain_len = numel(ff{1});
% for idx = 1:chain_len
%     FC{idx} = setItag(['K', suffixes(idx, :)], 'op', F);
% end

result = NRG_SL(nrgdata, H0, A0, Lambda, i * ff{1}, F, Z, 'Nkeep', Nkeep, 'dff', i * dff{1});

plotE(result, 'Emax', 10)
