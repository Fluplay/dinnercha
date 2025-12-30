function [Ek,Ek_patch] = get_Square(t,a,grid)

	%-2tcos(ka)
	Ek_patch = [];

	Ek = zeros(grid*2,grid*2);

	kx = linspace(0,pi/a,grid);
	ky = linspace(0,pi/a,grid);

	for x = (1:numel(kx))

		for y = (1:numel(ky))

			Ek_patch(x,y) = -2*t*cos(kx(x)) - 2*t*cos(ky(y));

		end

	end


	Ek(grid+1:2*grid, grid+1:2*grid) = Ek_patch; %(0,0) => (pi,pi)
	Ek(grid:-1:1, grid+1:2*grid) = Ek_patch; %(0,0) => (-pi,pi)
	Ek(grid:-1:1, grid:-1:1) = Ek_patch; %(0,0) => (-pi,-pi)
	Ek(grid+1:2*grid,grid:-1:1) = Ek_patch; %(0,0) => (pi,-pi)




end
