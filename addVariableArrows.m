function addVariableArrows(gx, varStruct)
% 在 geoaxes 周围添加整洁的箭头变量注释（不带方向符号）
% 变量名后显示贡献百分比，箭头粗细随贡献度变化

    pos = gx.Position; % [x y w h]
    xOffset = 0.05;
    yOffset = 0.05;
    minWidth = 1.0/2;
    maxWidth = 5.0/2;

    minWidth2 = 6.0/2;
    maxWidth2 = 16.0/2;
    % 默认颜色设置（可自定义）
    % dirColor = struct('top', 'b', 'left', 'k', 'bottom', 'g', 'right', [0.5 0 0.5]);
dirColor = struct( ...
    'top',    [230, 76, 60]/255,     ... % 红色
    'left',   [52, 152, 219]/255,    ... % 蓝色
    'bottom', [243, 156, 18]/255,    ... % 橙色
    'right',  [39, 174, 96]/255);        % 绿色


    function addArrow(fromX, fromY, toX, toY, txt, color, weight)
        lineWidth = minWidth + (maxWidth - minWidth) * weight;
        headSize = minWidth2 + (maxWidth2 - minWidth2) * weight;
        if all([fromX, toX, fromY, toY] >= 0) && all([fromX, toX, fromY, toY] <= 1)
            annotation('textarrow', [fromX, toX], [fromY, toY], ...
                'String', txt, ...
                'FontSize', 18, ...
                'Color', color, ...
                'LineWidth', lineWidth, ...
                'HeadWidth', headSize,...
                'FontWeight', 'normal', ...
                'HeadStyle', 'vback2');
        else
            warning('跳过越界注释: "%s"', txt);
        end
    end

    directions = {'top', 'left', 'bottom', 'right'};
    for d = 1:numel(directions)
        dir = directions{d};
        if isfield(varStruct, dir) && isfield(varStruct, [dir '_contr'])
            vars = varStruct.(dir);
            contr = varStruct.([dir '_contr']);
            n = numel(vars);
            max_val = max(contr);
            for i = 1:n
                norm_w = contr(i) ;%/ max_val;
                txt = sprintf('%s \n(%.1f%%)', vars{i}, contr(i));
                switch dir
                    case 'top'
                        x = pos(1) + (i - 0.5) * pos(3)/n;
                        fromY = min(pos(2) + pos(4) + yOffset, 0.99);
                        toY = pos(2) + pos(4);
                        addArrow(x, fromY, x, toY, txt, dirColor.(dir), norm_w);
                    case 'bottom'
                        x = pos(1) + (i - 0.5) * pos(3)/n;
                        fromY = max(pos(2) - yOffset, 0.01);
                        toY = pos(2);
                        addArrow(x, fromY, x, toY, txt, dirColor.(dir), norm_w);
                    case 'left'
                        y = pos(2) + (i - 0.5) * pos(4)/n;
                        fromX = max(pos(1) - xOffset, 0.01);
                        toX = pos(1);
                        addArrow(fromX, y, toX, y, txt, dirColor.(dir), norm_w);
                    case 'right'
                        y = pos(2) + (i - 0.5) * pos(4)/n;
                        fromX = min(pos(1) + pos(3) + xOffset, 0.99);
                        toX = pos(1) + pos(3);
                        addArrow(fromX, y, toX, y, txt, dirColor.(dir), norm_w);
                end
            end
        end
    end
end
