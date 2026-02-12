function diffJplot_SU2(obj,U)
	addpath('../')
	figure; hold on;
	fnames = 'Nothing'

	if U>1.2
		
		fnames =  ['../B*U=',num2str(U),'*D=_j153719*']

	elseif U<1.3
		fnames =  ['../B*U=',num2str(U),'*D=_j153584*']
	end

	mats = dir(fnames);
	colors = hsv(length(mats));
	colors


	J_vals = zeros(length(mats), 1);
	for k = 1:length(mats)
   		strVal = extractBetween(mats(k).name, 'Jl=', '_');
		J_vals(k) = str2double(strVal{1});
	end	

	[sorted_J, sortIdx] = sort(J_vals);	

	step = 1;

	idx = 1
	mat_sorted = mats(sortIdx);

	
	lw = 2


	for i = 1:step:length(mats)
		mname = mat_sorted(i).name ;
		load(['../',mname]);
		
		lab = extractBetween(mname,'Jl=','_Js');
		if strcmp(obj,'SE')
			SE_tot = sum(SE_arr,2);


			plot(ocont,-imag(SE_tot(:,1,end)),'color',colors(idx,:),'DisplayName',['J=',lab{1}],'LineWidth',lw);

		elseif strcmp(obj,'RSE')

			SE_tot = sum(SE_arr,2);
			plot(ocont,real(SE_tot(:,1,end)),'color',colors(idx,:),'DisplayName',['J=',lab{1}],'LineWidth',lw);

		elseif strcmp(obj,'A')

			A_tot = sum(A_arr,2);
			colors(i,:)
			plot(ocont,A_tot(:,1,end),'color',colors(idx,:),'DisplayName',['J=',lab{1}],'LineWidth',lw);

		end

		idx = idx +1

	end
	title(['U = ',num2str(U)]);
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
