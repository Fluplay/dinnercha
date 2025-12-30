function [ozin,RhoV2in] = get_Bethe(ozin,D)

	RhoV2in = [];
	for its = (1:numel(ozin))
		if abs(ozin(its)) >= D

			RhoV2in(end+1) = 0.0;

		else

			RhoV2in(end+1) = (D/ (2*pi) ) * sqrt(1 - ((ozin(its)/D).^2) );


		end



	end
end
