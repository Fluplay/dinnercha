function ZplotU()

    figure; hold on;
    addpath('../')
    fnames =  '../B*D=_j153559*'


    mats = dir(fnames);
    colors = jet(length(mats));

    u_vals = zeros(length(mats), 1);
    for k = 1:length(mats)
        strVal = extractBetween(mats(k).name, 'U=', '_');
        u_vals(k) = str2double(strVal{1});
    end

    [sorted_U, sortIdx] = sort(u_vals);

    mat_sorted = mats(sortIdx);
    step = 1;
    idx = 1

	Zlst =[];
	for i = 1:step:length(mats)-1
	    disp(i)
        mname = mat_sorted(i).name ;
        load(['../',mname]);

		zeroplus = size(ocont,1)/2+1;

		SE_tot = real(sum(SE_arr,2));
		diffRS = diff(SE_tot) ./ diff(ocont);

		Zlst(end+1) = 1/(1- diffRS(zeroplus));
		
	end

	scatter(sorted_U(1:end-1),Zlst,'filled')


end

