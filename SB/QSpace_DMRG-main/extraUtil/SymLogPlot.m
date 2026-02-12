function SymLogPlot(X, Y, color, varargin)
    
    % <Description>
    %
    % Symmetric Log Plot
    %
    % <Input>
    % X : [numeric vector / cell vector of numeric vectors]
    %       if X is a numeric vector, 
    %           X is a vector of x-coordinates of data points
    %       if X is a cell vector of numeric vectors,
    %           each cell elements is x-coordinates of data points 
    %           specified by color and legend
    %
    % Y : [numeric vector / cell vector of numeric vectors]
    %       if Y is a numeric vector, 
    %           Y is a vector of y-coordinates of data points
    %       if Y is a cell vector of numeric vectors,
    %           each cell elements is y-coordinates of data points 
    %           specified by color and legend
    %
    % color : [char / char array]
    %       if color is a char, 
    %           color is the color of the scatter plot
    %       if color is a char array,
    %           each element is color corresponding to each cell element of X and Y
    %
    % (The sizes of cell arrays X, Y and the lengths of vectors color, legend
    % must all be identical)
    %
    % <Option>
    % 'XScale',... : [char] x-axis scale. 'linear', 'log', or 'symlog' (Default: 'symlog')
    % 'YScale',... : [char] y-axis scale. 'linear', 'log', or 'symlog' (Default: 'linear')
    % 'Xmin',... : [numeric] logarithm of the minimum absolute value of the x-coordinate
    %           only used when 'XScale' is 'symlog' (Default: -10)
    % 'Ymin',... : [numeric] logarithm of the minimum absolute value of the y-coordinate
    %           only used when 'YScale' is 'symlog' (Default: -10)
    % 'style',... : [char / cell array of char] line style of marker of respective
    %                               data sets   (Default: '.')
    % 'LineWidth',... : [numeric] line width to be used if lines are plotted
    %                       (Default: 1)
    % 'MarkerSize',... : [numeric] marker size to be used if markers are plotted
    %                       (Default: 10)
    % 'XcontAtzero' : if used, the data is shown continuously near x = 0 (Default: not used)
    %                   (Only works when XScale is symlog and when the X data contains both positive and negative values)
    % 'YcontAtzero' : if used, the data is shown continuously near y = 0 (Default: not used)
    %                   (Only works when YScale is symlog and when the Y data contains both positive and negative values)
    % 'xzerowidth',... : [numeric] the ratio of distance between minimum absolute
    %                       value Xticks with respect to total x-axis length 
    %                       (Default: 0.08) (0<Xzerowidth<0.9)
    % 'yzerowidth',... : [numeric] the ratio of distance between minimum absolute
    %                       value Yticks with respect to total x-axis length
    %                       (Default: 0.08) (0<Yzerowidth<0.9)
    % 'legends',... : [char / char array] 
    %       if legend is a char,
    %           legend is the legend of the scatter plot
    %       if legend is a char array,
    %           each element is legend corresponding to each cell element of X and Y
    %
    % <Output>
    % Plot of the X and Y data

    %% checking inputs
    if isequal(class(X),'double')
        X = {X};
        lenX = 1;
    elseif isequal(class(X),'cell')
        lenX = numel(X);
    else
        error('ERR: X must be either a numeric vector or a cell vector of numeric vectors');
    end

    if isequal(class(Y),'double')
        Y = {Y};
        lenY = 1;
    elseif isequal(class(Y),'cell')
        lenY = numel(Y);
    else
        error('ERR: Y must be either a numeric vector or a cell vector of numeric vectors');
    end
    
    if isequal(class(color),'char')
        color = {color};
        lenC = 1;
    elseif isequal(class(color),'cell')
        lenC = numel(color);
    else
        error('ERR: color must be a char or a cell vector of char');
    end

    if ~isequal(lenX, lenY, lenC)
        error('ERR: The number of data types, colors, and legends must be consistent');
    end


    %% Parsing options

    % Default values of options
    XScale = 'symlog';
    YScale = 'linear';
    Xmin = -10;
    Ymin = -10;
    style = cell(numel(X),1);
    for it = (1:numel(X))
        style{it} = '.';
    end
    linewidth = 1;
    MarkerSize = 10;
    Xzerowidth = 0.08;
    Yzerowidth = 0.08;
    XcontAtzero = false;
    YcontAtzero = false;
    PlotLegend = false;

    while ~isempty(varargin)
        switch varargin{1}
            case 'XScale'
                if ismember(varargin{2}, {'linear', 'log', 'symlog'})
                    XScale = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''XScale'' must be ''linear'', ''log'', or ''symlog''');
                end

            case 'YScale'
                if ismember(varargin{2}, {'linear', 'log', 'symlog'})
                    YScale = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''XScale'' must be ''linear'', ''log'', or ''symlog''');
                end

            case 'Xmin'
                if isnumeric(varargin{2})
                    if varargin{2} ~= 0
                        XScale = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: ''Xmin'' cannot be zero');
                    end
                else
                    error('ERR: ''Xmin'' must be a number');
                end

            case 'Ymin'
                if isnumeric(varargin{2})
                    if varargin{2} ~= 0
                        YScale = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: ''Ymin'' cannot be zero');
                    end
                else
                    error('ERR: ''Ymin'' must be a number');
                end

            case 'style'
                if isequal(class(varargin{2}),'char')
                    style = {varargin{2}};
                    varargin(1:2) = [];
                elseif isequal(class(varargin{2}),'cell')
                    style = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: unknown input type for option ''style''');
                end

            case 'linewidth'
                if isnumeric(varargin{2})
                    if varargin{2} > 0
                        linewidth = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: linewidth must be a positive real number');
                    end
                else
                    error('ERR: unknown input type for option ''linewidth''');
                end

            case 'LineWidth'
                if isnumeric(varargin{2})
                    if varargin{2} > 0
                        linewidth = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: linewidth must be a positive real number');
                    end
                else
                    error('ERR: unknown input type for option ''linewidth''');
                end

            case 'MarkerSize'
                if isnumeric(varargin{2})
                    if varargin{2} > 0
                        MarkerSize = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: MarkerSize must be a positive real number');
                    end
                else
                    error('ERR: unknown input type for option ''MarkerSize''');
                end

            case 'xzerowidth'
                if isnumeric(varargin{2})
                    if varargin{2} >= 0 && varargin{2} <0.9
                        Xzerowidth = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: xzerowidth must be a between 0 and 0.9');
                    end
                else
                    error('ERR: xzerowidth must be a number');
                end

            case 'yzerowidth'
                if isnumeric(varargin{2})
                    if varargin{2} >= 0 && varargin{2} <0.9
                        Yzerowidth = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: yzerowidth must be a nonnegative real number smaller than 0.9');
                    end
                else
                    error('ERR: yzerowidth must be a number');
                end

            case 'legends'
                PlotLegend = true;
                if isequal(class(varargin{2}),'char')
                    legends = {varargin{2}};
                    varargin(1:2) = [];
                    lenL = 1;
                elseif isequal(class(varargin{2}),'cell')
                    legends = varargin{2};
                    varargin(1:2) = [];
                    lenL = numel(legends);
                else
                    error('ERR: legend must be a char or a cell vector of char');
                end

                if lenL ~= lenX
                    error('ERR: The number of legends must be consistent with the number of data types, colors, and legends');
                end

            case 'XcontAtzero'
                XcontAtzero = true;
                varargin(1) = [];
                Xzerowidth = 0;

            case 'YcontAtzero'
                YcontAtzero = true;
                varargin(1) = [];
                Yzerowidth = 0;

            otherwise
                if ischar(varargin{1})
                    error(['ERR: Unknown option ',varargin{1}]);
                else
                    error('ERR: Unknown option');
                end
        end
    end

    if XcontAtzero && ~isequal(XScale,'symlog')
        error('''XcontAtzero'' can only be used when XScale is symlog')
    end

    if YcontAtzero && ~isequal(YScale,'symlog')
        error('''YcontAtzero'' can only be used when XScale is symlog')
    end

    %% Convert data for plotting
    if isequal(XScale,'symlog')

        sign_X = '';
        minExp_arr = zeros(numel(X), 1);
        maxExp_arr = zeros(numel(X), 1);

        for it = (1:numel(X))

            X_temp = X{it};

            if ~all(X_temp~=0)
                error('ERR: X must not contain zero in symlog scale');
            end

            if all(X_temp>0)
                minExp_arr(it) = floor(log(min(X_temp))/log(10));
                maxExp_arr(it) = ceil(log(max(X_temp))/log(10));
                if ~isequal(sign_X,'undefined')
                    if isequal(sign_X,'negative')
                        sign_X = 'undefined';
                    else
                        sign_X = 'positive';
                    end
                end

            elseif all(X_temp<0)
                minExp_arr(it) = floor(log(-max(X_temp))/log(10));
                maxExp_arr(it) = ceil(log(-min(X_temp))/log(10));
                if ~isequal(sign_X,'undefined')
                    if isequal(sign_X,'positive')
                        sign_X = 'undefined';
                    else
                        sign_X = 'negative';
                    end
                end

            else
                pos_X = X_temp(X_temp>0);
                neg_X = X_temp(X_temp<0);
                minExp_arr(it) = min( floor(log(min(pos_X))/log(10)), floor(log(-max(neg_X))/log(10)) );
                maxExp_arr(it) = max( ceil(log(max(pos_X))/log(10)), ceil(log(-min(neg_X))/log(10)) );
                sign_X = 'undefined';

            end
        end

        minExp = min(minExp_arr);
        maxExp = max(maxExp_arr);

        if maxExp - minExp > 5
            minTick = 5*floor(minExp/5);
            maxTick = 5*floor(maxExp/5);
            tickInterval = 5*ceil((maxTick-minTick)/25);
        else
            minTick = minExp;
            maxTick = maxExp;
            tickInterval = 1;
        end

        % Normalize x-coordinates
        for Idx = (1:numel(X))

            X_temp = X{Idx};
            X_norm = zeros(numel(X_temp), 1);
            for it = (1:numel(X_temp))
                X_norm(it) = Xzerowidth + (1-Xzerowidth-0.05)*(log(abs(X_temp(it)))/log(10) - minTick)/(maxExp - minTick);
                X_norm(it) = sign(X_temp(it))*X_norm(it);
            end
            X{Idx} = X_norm;
        end

        % generate xtick positions and their labels
        temp = (minTick : tickInterval : maxTick);
        XtickPositions = [];
        for it = (1:numel(temp))
            XtickPositions = [XtickPositions, ...
                                Xzerowidth + (1-Xzerowidth-0.05)*(temp(it) - minTick)/(maxExp - minTick)];
        end
        
        XtickLabels = cell(0,0);
        switch sign_X
            case 'undefined'
                if XcontAtzero
                    XtickPositions = [-flip(XtickPositions(2:end)),XtickPositions];
                    
                    temp = flip(temp);
                    for it = (1:numel(temp)-1)
                        XtickLabels = cat(2,XtickLabels,'$-10^{'+string(temp(it))+'}$');
                    end
                    XtickLabels = cat(2,XtickLabels,'$\pm10^{'+string(temp(end))+'}$');
                    temp = flip(temp);
                    for it = (2:numel(temp))
                        XtickLabels = cat(2,XtickLabels,'$10^{'+string(temp(it))+'}$');
                    end
                else
                    XtickPositions = [-flip(XtickPositions), 0, XtickPositions];

                    temp = flip(temp);
                    for it = (1:numel(temp))
                        XtickLabels = cat(2,XtickLabels,'$-10^{'+string(temp(it))+'}$');
                    end
                    XtickLabels = cat(2,XtickLabels,'0');
                    temp = flip(temp);
                    for it = (1:numel(temp))
                        XtickLabels = cat(2,XtickLabels,'$10^{'+string(temp(it))+'}$');
                    end
                end

                

            case 'positive'

                if XcontAtzero
                    error('''XcontAtzero'' can only be used when the X data contains both positive and negative values');
                end
                XtickPositions = [0, XtickPositions];

                for it = (1:numel(temp))
                    XtickLabels = cat(2,XtickLabels,'$10^{'+string(temp(it))+'}$');
                end
                XtickLabels = cat(2,'0',XtickLabels);

            case 'negative'

                if XcontAtzero
                    error('''XcontAtzero'' can only be used when the X data contains both positive and negative values');
                end
                XtickPositions = [-flip(XtickPositions), 0];

                temp = flip(temp);
                for it = (1:numel(temp))
                    XtickLabels = cat(2,XtickLabels,'$-10^{'+string(temp(it))+'}$');
                end
                XtickLabels = cat(2,XtickLabels,'0');
        end

    end


    if isequal(YScale,'symlog')

        sign_Y = '';
        minExp_arr = zeros(numel(X), 1);
        maxExp_arr = zeros(numel(X), 1);

        for it = (1:numel(Y))

            Y_temp = Y{it};

            if ~all(Y_temp~=0)
                error('ERR: Y must not contain zero in symlog scale');
            end

            

            if all(Y_temp>0)
                minExp_arr(it) = floor(log(min(Y_temp))/log(10));
                maxExp_arr(it) = ceil(log(max(Y_temp))/log(10));
                if ~isequal(sign_Y,'undefined')
                    if isequal(sign_Y,'negative')
                        sign_Y = 'undefined';
                    else
                        sign_Y = 'positive';
                    end
                end

            elseif all(Y_temp<0)
                minExp_arr(it) = floor(log(-max(Y_temp))/log(10));
                maxExp_arr(it) = ceil(log(-min(Y_temp))/log(10));
                if ~isequal(sign_Y,'undefined')
                    if isequal(sign_Y,'positive')
                        sign_Y = 'undefined';
                    else
                        sign_Y = 'negative';
                    end
                end

            else
                pos_Y = Y_temp(Y_temp>0);
                neg_Y = Y_temp(Y_temp<0);
                minExp_arr(it) = min( floor(log(min(pos_Y))/log(10)), floor(log(-max(neg_Y))/log(10)) );
                maxExp_arr(it) = max( ceil(log(max(pos_Y))/log(10)), ceil(log(-min(neg_Y))/log(10)) );
                sign_Y = 'undefined';
            end
        end

        minExp = min(minExp_arr);
        maxExp = max(maxExp_arr);

        if maxExp - minExp > 5
            minTick = 5*floor(minExp/5);
            maxTick = 5*floor(maxExp/5);
            tickInterval = 5*ceil((maxTick-minTick)/25);
        else
            minTick = minExp;
            maxTick = maxExp;
            tickInterval = 1;
        end

        % Normalize y-coordinates
        for Idx = (1:numel(Y))

            Y_temp = Y{Idx};
            Y_norm = zeros(numel(Y_temp), 1);
            for it = (1:numel(Y_temp))
                Y_norm(it) = Yzerowidth + (1-Yzerowidth-0.05)*(log(abs(Y_temp(it)))/log(10) - minTick)/(maxExp - minTick);
                Y_norm(it) = sign(Y_temp(it))*Y_norm(it);
            end
            Y{Idx} = Y_norm;
        end

        % generate ytick positions and their labels
        temp = (minTick : tickInterval : maxTick);
        YtickPositions = [];
        for it = (1:numel(temp))
            YtickPositions = [YtickPositions, ...
                                Yzerowidth + (1-Yzerowidth-0.05)*(temp(it) - minTick)/(maxExp - minTick)];
        end

        YtickLabels = cell(0,0);
        switch sign_Y
            case 'undefined'
                if YcontAtzero
                    YtickPositions = [-flip(YtickPositions(2:end)),YtickPositions];

                    temp = flip(temp);
                    for it = (1:numel(temp)-1)
                        YtickLabels = cat(2,YtickLabels,'$-10^{'+string(temp(it))+'}$');
                    end
                    YtickLabels = cat(2,YtickLabels,'$\pm10^{'+string(temp(end))+'}$');
                    temp = flip(temp);
                    for it = (2:numel(temp))
                        YtickLabels = cat(2,YtickLabels,'$10^{'+string(temp(it))+'}$');
                    end
                else
                    YtickPositions = [-flip(YtickPositions), 0, YtickPositions];

                    temp = flip(temp);
                    for it = (1:numel(temp))
                        YtickLabels = cat(2,YtickLabels,'$-10^{'+string(temp(it))+'}$');
                    end
                    YtickLabels = cat(2,YtickLabels,'0');
                    temp = flip(temp);
                    for it = (1:numel(temp))
                        YtickLabels = cat(2,YtickLabels,'$10^{'+string(temp(it))+'}$');
                    end
                end

            case 'positive'

                if YcontAtzero
                    error('''YcontAtzero'' can only be used when the Y data contains both positive and negative values');
                end
                YtickPositions = [0, YtickPositions];

                for it = (1:numel(temp))
                    YtickLabels = cat(2,YtickLabels,'$10^{'+string(temp(it))+'}$');
                end
                YtickLabels = cat(2,'0',YtickLabels);

            case 'negative'

                if YcontAtzero
                    error('''YcontAtzero'' can only be used when the Y data contains both positive and negative values');
                end
                YtickPositions = [-flip(YtickPositions), 0];

                temp = flip(temp);
                for it = (1:numel(temp))
                    YtickLabels = cat(2,YtickLabels,'$-10^{'+string(temp(it))+'}$');
                end
                YtickLabels = cat(2,YtickLabels,'0');
        end        
        
        [YtickPositions, ids] = sort(YtickPositions);
        YtickLabels = YtickLabels(ids);
    end


%% Plot
hold on;
    if isequal(XScale,'symlog')
        if isequal(YScale,'symlog')

            for it = (1:numel(X))
                if ismember(style{it},{'-', '--', ':', '-.'})
                    plot(X{it},Y{it},style{it},'LineWidth',linewidth,'color',color{it});
                else
                    plot(X{it},Y{it},style{it},'MarkerSize',MarkerSize,'color',color{it});
                end
            end
            set(gca,'XScale','linear','YScale','linear');
            XAxisProperties = get(gca, 'XAxis');
            XAxisProperties.FontSize = 15;
            XAxisProperties.TickLabelInterpreter = 'latex';
            YAxisProperties = get(gca, 'YAxis');
            YAxisProperties.FontSize = 15;
            YAxisProperties.TickLabelInterpreter = 'latex';
            % For conventional uses, 'tex'

            if isequal(sign_X,'negative')
                xlim([-1,0.1]);
            elseif isequal(sign_X,'positive')
                xlim([-0.1,1]);
            else
                xlim([-1,1]);
            end

            if isequal(sign_Y,'negative')
                ylim([-1,0.1]);
            elseif isequal(sign_Y,'positive')
                ylim([-0.1,1]);
            else
                ylim([-1,1]);
            end

            xticks(XtickPositions);
            xticklabels(XtickLabels);
            yticks(YtickPositions);
            yticklabels(YtickLabels);
            if PlotLegend
                legend(legends,'Location','northeast','FontSize',15);
            end

        else
            
            for it = (1:numel(X))
                if ismember(style{it},{'-', '--', ':', '-.'})
                    plot(X{it},Y{it},style{it},'LineWidth',linewidth,'color',color{it});
                else
                    plot(X{it},Y{it},style{it},'MarkerSize',MarkerSize,'color',color{it});
                end
            end
            set(gca,'XScale','linear','YScale',YScale);
            XAxisProperties = get(gca, 'XAxis');
            XAxisProperties.FontSize = 15;
            XAxisProperties.TickLabelInterpreter = 'latex';
            YAxisProperties = get(gca, 'YAxis');
            YAxisProperties.FontSize = 15;
            YAxisProperties.TickLabelInterpreter = 'latex';
            % For conventional uses, 'tex'
            xlim([-1,1]);
            xticks(XtickPositions);
            xticklabels(XtickLabels);
            if PlotLegend
                legend(legends,'Location','northeast','FontSize',15);
            end

        end

    else
        if isequal(YScale,'symlog')
            for it = (1:numel(X))
                if ismember(style{it},{'-', '--', ':', '-.'})
                    plot(X{it},Y{it},style{it},'LineWidth',linewidth,'color',color{it});
                else
                    plot(X{it},Y{it},style{it},'MarkerSize',MarkerSize,'color',color{it});
                end
            end
            set(gca,'XScale',XScale,'YScale','linear');
            XAxisProperties = get(gca, 'XAxis');
            XAxisProperties.FontSize = 15;
            XAxisProperties.TickLabelInterpreter = 'latex';
            YAxisProperties = get(gca, 'YAxis');
            YAxisProperties.FontSize = 15;
            YAxisProperties.TickLabelInterpreter = 'latex'; 
            % For conventional uses, 'tex'
            if isequal(sign_Y,'negative')
                ylim([-1,0.1]);
            elseif isequal(sign_Y,'positive')
                ylim([-0.1,1]);
            else
                ylim([-1,1]);
            end
            yticks(YtickPositions);
            yticklabels(YtickLabels);
            if PlotLegend
                legend(legends,'Location','northeast','FontSize',15);
            end

        else
            for it = (1:numel(X))
                if ismember(style{it},{'-', '--', ':', '-.'})
                    plot(X{it},Y{it},style{it},'LineWidth',linewidth,'color',color{it});
                else
                    plot(X{it},Y{it},style{it},'MarkerSize',MarkerSize,'color',color{it});
                end
            end
            set(gca,'XScale',XScale,'YScale',YScale);
            if PlotLegend
                legend(legends,'Location','northeast','FontSize',15);
            end
        end

    end

    
end
