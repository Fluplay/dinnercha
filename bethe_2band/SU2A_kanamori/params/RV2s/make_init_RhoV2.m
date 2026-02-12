function make_init_RhoV2(fname,fo_or_ba)

	prev_params = load(fname);

	RV2arr = prev_params.RhoV2in_arr;

	RhoV2prev = RV2arr(:,:,end);
	
	U = prev_params.U;

	save_name = ['U',num2str(U),'_',fo_or_ba,'.mat'];

	save(save_name,'RhoV2prev');

	disp([save_name,' is saved '])

end
