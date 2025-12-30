function Ek_patch = get_Cubic(t,a,grid)


    %-2tcos(ka)
    Ek_patch = [];

    Ek = zeros(grid*2,grid*2);

    kx = linspace(0,pi/a,grid);
    ky = linspace(0,pi/a,grid);
    kz = linspace(0,pi/a,grid)

    for x = (1:numel(kx))

        for y = (1:numel(ky))

            for z = (1:numel(kz))

                Ek_patch(x,y,z) = -2*t*cos(kx(x)) - 2*t*cos(ky(y)) -2*t*cos(kz(z));

            end


        end

    end




end
