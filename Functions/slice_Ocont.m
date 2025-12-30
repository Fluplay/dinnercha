function oconts = slice_Ocont(ocont,sn)
	
	oconts = {}

	piece = floor(numel(ocont)/sn)

	for i = (1:sn)
		if i == sn

			oconts{i} = ocont(1 + piece*(i-1):end)
		else

			oconts{i} = ocont(1 + piece*(i-1) : piece*i )
		end
	end

end
		

