function A = getRandDataSB(A)
    for i = 1:numel(A.data)
        A.data{i} = rand(size(A.data{i}));
    end
end