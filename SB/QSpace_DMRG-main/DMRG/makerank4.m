function rank4op = makerank4(operator)
    if rank(operator) == 2
        rank4op = appendSingletons(operator, '  +-'); % lower, upper, left, right
    elseif rank(operator) == 3
        if ~isempty(operator.info.itags{3}) && operator.info.itags{3}(end)=='*' % out
            rank4op = permute(appendSingletons(operator, '+'), [1 2 4 3]); % lower, upper, left, right
        else % in
            rank4op = appendSingletons(operator, '-'); % lower, upper, left, right
        end
    elseif rank(operator) == 4
        rank4op=operator;
    else
        error('please check the rank of operator in makerank4')
    end
end
