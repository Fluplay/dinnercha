function plot_convergency(path,obj)
	
	last_N = 5;

	load(path);


	P = 'nothing';

	if obj == 'A'
	
		P = A_arr;

	elseif obj == 'R'

		P = RhoV2in_arr;

	elseif obj == 'SE'

		P = SE_arr;

	else
		error('wrong obj')

	end

	if length(size(P)) ~= 3
		error('Wrong size check DMFT step')
	end

	pobj = sum(P,2);

	figure; hold on;


	if obj == 'SE'
		pobj = -imag(pobj)
	end

	Nstep = size(P,3)

	colors = jet(length(( Nstep:-1:(Nstep-last_N+1))))

	if Nstep < last_N
		last_N = Nstep
	end
	
	color_idx = 1
	for i = Nstep:-1:(Nstep-last_N+1)
		plot(ocont,pobj(:,1,i),'color',colors(color_idx,:),'DisplayName',['Step',num2str(i)])
		color_idx = color_idx +1
	end
	
	if obj =='SE'

        set(gca,'XScale','log')
        set(gca,'YScale','log')
        xlim([1e-6,5])

    else

        xlim([-5,5])

    end
    legend('show')

end
