function [xi,Rsquare] = exponentRegression(x,y)
    if min(y)<= 0
        % error('y must be positive');
        xi = NaN;
        Rsquare = NaN;
        return
    end

    mdl = fitlm(x,log(y));
    slope = mdl.Coefficients.Estimate(2);
    Rsquare = mdl.Rsquared.Ordinary;
    xi = -1/slope;
end