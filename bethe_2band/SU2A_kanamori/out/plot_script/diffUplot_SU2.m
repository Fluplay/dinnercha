function diffUplot_SU2(obj)
	addpath('../')
	figure; hold on;
	fnames =  '../B*D=_j153559*'


	mats = dir(fnames);
	colors = lines(length(mats));

	u_vals = zeros(length(mats), 1);
	for k = 1:length(mats)
   		strVal = extractBetween(mats(k).name, 'U=', '_');
		u_vals(k) = str2double(strVal{1});
	end	

	[sorted_U, sortIdx] = sort(u_vals);	

	mat_sorted = mats(sortIdx);
	step = 1;
	idx = 1
	lw = 2
	for i = 1:step:length(mats)-1
		disp(i)
		mname = mat_sorted(i).name ;
		load(['../',mname]);
		
		lab = extractBetween(mname,'jljs_','_Jl');
		%color_idx = str2num(lab{1}(3:end))*5 + 1
		if strcmp(obj,'SE')
			SE_tot = sum(SE_arr,2);


			plot(ocont,-imag(SE_tot(:,1,end)),'color',colors(idx,:),'DisplayName',lab{1},'LineWidth',lw);



		elseif strcmp(obj,'RSE')

			SE_tot = sum(SE_arr,2);

			plot(ocont,real(SE_tot(:,1,end)),'color',colors(idx,:),'DisplayName',lab{1},'LineWidth',lw);


		elseif strcmp(obj,'A')

			A_tot = sum(A_arr,2);

			plot(ocont,A_tot(:,1,end),'color',colors(idx,:),'DisplayName',lab{1},'LineWidth',lw);

		end

		idx = idx +1 
	end
	if strcmp(obj,'SE')

		set(gca,'XScale','log')
		set(gca,'YScale','log')


	elseif strcmp(obj,'RSE')
		xlim([-10,10])


	elseif strcmp(obj,'A')

		xlim([-5,5])

	end
	legend('show')



end
