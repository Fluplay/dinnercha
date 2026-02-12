function Mresult = dividePhysicalBond(M,div1,div2)
% div1, div2: index of quantum number to divide
    Nsite = numel(M);
    Mresult = QSpace(Nsite,1);
    for itN = 1:Nsite
        Ms = M(itN);
        Ms = appendSingletons(Ms,'+');
        Ms.Q{4} = Ms.Q{3};
        Ms.Q{3}(:,div2) = zeros(size(Ms.Q{3}(:,div2)));
        Ms.Q{4}(:,div1) = zeros(size(Ms.Q{4}(:,div1)));
        for k = 1:size(Ms.Q{4},2)
            if k ~= div1 && k ~= div2
                Ms.Q{4}(:,k) = zeros(size(Ms.Q{4}(:,k)));
            end
        end
        for i = 1:numel(Ms.data)
            Ms.info.cgr(i,div2).qset =  Ms.info.cgr(i,div2).qset([1 2 4 3]);
            Ms.info.cgr(i,div2).size =  Ms.info.cgr(i,div2).size([1 2 4 3]);
        end
        Mresult(itN) = Ms;
    end

end