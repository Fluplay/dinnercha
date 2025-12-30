function G = get_GF0(ocont,Ek,eta)

    G = [];
    for ome_idx = (1:numel(ocont))
        disp("Processing to get GF -> " + ome_idx/numel(ocont) *100)
        Gtmp = 0.0 + 0.0j;

        for i = (1:numel(Ek))

            Gtmp = Gtmp + (1 / (ocont(ome_idx) - Ek(i) + eta*1j));

        end
        G(end + 1) = Gtmp ;
    end
    G = G/numel(Ek);

	G = zeros(numel(ocont),1);
	for om = (1:numel(ocont))
		G(om) = ocont(om)+eta*1j;
		for k = (1:numel(Ek))
			G(om) = G(om)+1/(G(om)-Ek(k));
		end
	end

end

