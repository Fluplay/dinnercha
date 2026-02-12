
function ZplotJ(U)

    addpath('../')
    figure; hold on;
    fnames = 'Nothing'

    if U>1.2

        fnames =  ['../B*U=',num2str(U),'*D=_j153719*']

    elseif U<1.3
        fnames =  ['../B*U=',num2str(U),'*D=_j153584*']
    end

    mats = dir(fnames);
    colors = turbo(length(mats));
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



    Zlst =[];
    for i = 1:step:length(mats)
        disp(i)
        mname = mat_sorted(i).name ;
        load(['../',mname]);

        zeroplus = size(ocont,1)/2+2;

        SE_tot = real(sum(SE_arr,2));
        diffRS = diff(SE_tot) ./ diff(ocont);

        Zlst(end+1) = 1/(1- diffRS(zeroplus));

    end

    scatter(sorted_J,Zlst,'filled')
	title(['U = ',num2str(U)]);
end
