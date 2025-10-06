function fig=plotStackedBarWithTruePies(years, stackData, pieData, barColors, pieColors)

    xPos = 1:length(years);
    fig=figure('Color','w');
    set(gcf,'Position',[ 563   278   784   497])
    ax = axes;
    hold on;

    % Plot stacked bars
    b = bar(xPos, stackData', 'stacked');
    for i = 1:length(b)
        b(i).FaceColor = barColors(i,:);
    end

    % Axis settings
    xticks(xPos);
    xticklabels(string(years));
    ylabel('Contribution (%)');
    set(gca, 'FontSize', 12,'FontWeight','bold');
    ylim([0, max(sum(stackData,1)) * 1.3]);
    yticks([0:20:100]);
    % box on
    % set(gca,'fontsize',12)

    % Legend
    % legend(b, strcat("Cat ", string(1:size(stackData,1))), ...
    %        'Location','eastoutside', 'Box','off');
    le=legend(b,'PRE', 'EDD', 'GDD','NT', 'IRR', 'NFER', 'DRA', 'CC','BD', 'CLAY', 'SAND', 'SILT', 'PH','Location','east', 'Box','off');
    le.Position=[0.8122    0.1449    0.1298    0.5171];

    % Get axes and figure position info
    axPos = get(ax, 'Position');
    xLimits = xlim(ax);
    yLimits = ylim(ax);

    % Place pie charts
    for i = 1:length(xPos)
        xData = xPos(i);
        yData = sum(stackData(:,i)) + max(sum(stackData,1)) * 0.1;

        % --- Convert (xData, yData) from data coords → normalized figure units
        xNorm = (xData - xLimits(1)) / (xLimits(2) - xLimits(1));
        yNorm = (yData - yLimits(1)) / (yLimits(2) - yLimits(1));
        xFig = axPos(1) + xNorm * axPos(3);
        yFig = axPos(2) + yNorm * axPos(4);

        % Create tiny axes for pie chart at the correct normalized position
        pieSize = 0.1;
        axPie = axes('Position', [xFig - pieSize/2, yFig - pieSize/2, pieSize, pieSize]);
        % labels={'C','M','S'};
        
        p=pie(axPie, pieData(i,:));
        colormap(axPie, pieColors);
        set(axPie, 'Visible', 'off');        
        
        

    end

    % title('Stacked Bar with Always-Round Pie Charts');
end
