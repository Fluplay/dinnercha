function diffaLplot(obj)
	figure; hold on;
	fnames = '../*U=1.4*'
	mats = dir(fnames);
	colors = cool(length(mats));

    u_vals = zeros(length(mats), 1);
    for k = 1:length(mats)
        strVal = extractBetween(mats(k).name, 'alphaL=', '_');
        aL_vals(k) = str2double(strVal{1});
    end

    [sorted_aL, sortIdx] = sort(aL_vals);

    mat_sorted = mats(sortIdx)
	for i = 1:length(mats)
		mname = mat_sorted(i).name ;
		load(['../',mname]);
		
		lab = extractBetween(mname,'nz=2_','_nd');
		color_idx = str2num(lab{1}(8:end))*10 + 1;
		if obj == 'SE'
			SE_tot = sum(SE_arr,2);

			plot(ocont,-imag(SE_tot(:,1,end)),'color',colors(color_idx,:),'DisplayName',lab{1});

		


		elseif obj == 'A'

			A_tot = sum(A_arr,2);

			plot(ocont,A_tot(:,1,end),'color',colors(color_idx,:),'DisplayName',lab{1});

		end

	end

	if obj == 'SE'

		set(gca,'XScale','log')
		set(gca,'YScale','log')
		xlim([1e-6,5])

	elseif obj == 'A'

		xlim([-5,5])

	end
	legend('show')



end
