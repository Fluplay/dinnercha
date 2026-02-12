function diffUplot(obj)
	figure; hold on;
	fnames =  '../*alphaL=1_nd=2*'
	mats = dir(fnames);
	colors = turbo(length(mats));

	u_vals = zeros(length(mats), 1);
	for k = 1:length(mats)
   		strVal = extractBetween(mats(k).name, 'U=', '_');
		u_vals(k) = str2double(strVal{1});
	end	

	[sorted_U, sortIdx] = sort(u_vals);	

	mat_sorted = mats(sortIdx)

	for i = 1:length(mats)
		mname = mat_sorted(i).name ;
		load(['../',mname]);
		
		lab = extractBetween(mname,'jljs_','_Jl');
		color_idx = str2num(lab{1}(3:end))*5 + 1;
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
