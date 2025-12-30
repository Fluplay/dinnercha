function [ocont,G] = glue_Ocont(Path)
	
	pieces = [Path,'GF0_*'];
	merged = [Path,'Merged.mat'];
	files = dir(pieces);
	ocont_cell = {};
	G_cell = {};

	disp("trying to merge sliced data => 1 to " + numel(files))
	for i = 1:numel(files)
		
		file_name = [Path,files(i).name];

		disp("processing  -> " + i + "/" + numel(files)) 

		load(file_name);

		ocont_cell{i,1} = ocont;

		G_cell{i,1} = reshape(G,[],1);
	end

	ocont_merged = cell2mat(ocont_cell);

	G_merged = cell2mat(G_cell);

	ocont = ocont_merged;

	[ocont,indecies] = sort(ocont_merged);

	G = G_merged(indecies);

	clear G_merged ocont_merged;

	save(merged)


end


		
