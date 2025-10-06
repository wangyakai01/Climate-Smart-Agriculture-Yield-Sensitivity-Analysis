%% ===============================================================
%       FIGURE 1: Climate + Management + Soil Effects
%       Analysis of Food Production Components (FEU / Energy / Protein)
% ===============================================================

clear; clc;

%% ===============================================================
%   Load data
loadAllVars('.\Data\')
%% ===============================================================
%   STEP 1. Compute Combined Food, Energy, and Protein Production
% ===============================================================

% ---- Nutritional coefficients for each crop ----
Energy_kJ = [1527,1590,1368,1377,1548,1418,1473];  % Energy content (kJ per unit yield)
Protein_g = [9.42,34.1,12.61,15.4,7.94,11.3,9.91]; % Protein content (g per unit yield)

% Weighting coefficients for Food Equivalent Unit (FEU)
C_H = 0.9 / 1527;   % Energy weighting
C_P = 0.1 / 9.42;   % Protein weighting
FEU_coef     = C_H .* Energy_kJ + C_P .* Protein_g;
Energy_coef  = C_H .* Energy_kJ;
Protein_coef = C_P .* Protein_g;

% Crop area dataset names
areanames = {'cornareadata1','soybeanareadata1','wheat_winterareadata1', ...
             'wheat_springareadata1','RICEareadata1','SORGHUMareadata1','BARLEYareadata1'};

% ---- Combine yield × area × conversion coefficient ----
Allfeudata = cat(3, ...
    cornyielddata      * 62.77 * FEU_coef(1)     .* cornareadata1, ...
    soybeanyielddata   * 67.25 * FEU_coef(2)     .* soybeanareadata1, ...
    wheat_winteryielddata * 67.25 * FEU_coef(3) .* wheat_winterareadata1, ...
    wheat_springyielddata * 67.25 * FEU_coef(4) .* wheat_springareadata1, ...
    RICEyielddata      * 1.121 * FEU_coef(5)     .* RICEareadata1, ...
    SORGHUMyielddata   * 62.77 * FEU_coef(6)     .* SORGHUMareadata1, ...
    BARLEYyielddata    * 53.8  * FEU_coef(7)     .* BARLEYareadata1);

Allenergydata = cat(3, ...
    cornyielddata * 62.77 * Energy_coef(1) .* cornareadata1, ...
    soybeanyielddata * 67.25 * Energy_coef(2) .* soybeanareadata1, ...
    wheat_winteryielddata * 67.25 * Energy_coef(3) .* wheat_winterareadata1, ...
    wheat_springyielddata * 67.25 * Energy_coef(4) .* wheat_springareadata1, ...
    RICEyielddata * 1.121 * Energy_coef(5) .* RICEareadata1, ...
    SORGHUMyielddata * 62.77 * Energy_coef(6) .* SORGHUMareadata1, ...
    BARLEYyielddata * 53.8 * Energy_coef(7) .* BARLEYareadata1);

Allproteindata = cat(3, ...
    cornyielddata * 62.77 * Protein_coef(1) .* cornareadata1, ...
    soybeanyielddata * 67.25 * Protein_coef(2) .* soybeanareadata1, ...
    wheat_winteryielddata * 67.25 * Protein_coef(3) .* wheat_winterareadata1, ...
    wheat_springyielddata * 67.25 * Protein_coef(4) .* wheat_springareadata1, ...
    RICEyielddata * 1.121 * Protein_coef(5) .* RICEareadata1, ...
    SORGHUMyielddata * 62.77 * Protein_coef(6) .* SORGHUMareadata1, ...
    BARLEYyielddata * 53.8 * Protein_coef(7) .* BARLEYareadata1);

Allareadata = cat(3, cornareadata1, soybeanareadata1, wheat_winterareadata1, ...
                     wheat_springareadata1, RICEareadata1, SORGHUMareadata1, BARLEYareadata1);

% ---- Sum across crops ----
Feu_sum     = nansum(Allfeudata, 3);
Energy_sum  = nansum(Allenergydata, 3);
Protein_sum = nansum(Allproteindata, 3);
AREA_sum    = nansum(Allareadata, 3);

% Avoid division by zero
AREA_sum(AREA_sum==0) = nan;
Feu_sum(Feu_sum==0)   = nan;
Energy_sum(Energy_sum==0) = nan;
Protein_sum(Protein_sum==0) = nan;

% ---- Mean per area ----
Feu_mean     = Feu_sum ./ AREA_sum;
Energy_mean  = Energy_sum ./ AREA_sum;
Protein_mean = Protein_sum ./ AREA_sum;


%% ===============================================================
%   STEP 2. Random Forest: Effects of Climate, Management, Soil
% ===============================================================

precall = metindexall(1).prec;
GDD_sum = metindexall(1).gdd;
EDD_sum = metindexall(1).edd;

Period_range = [1:10; 11:20; 21:30; 31:40];  % 4 decades (1980s–2010s)
predictorNames = {'PRE','EDD','GDD','NT','IRR','NFER','TILE','CC','bd','clay','sand','silt','ph'};

% Storage variables
Feu_p = nan(3109,4); Energy_p = nan(3109,4); Protein_p = nan(3109,4);
impall_feu = nan(4,13); impall_energy = nan(4,13); impall_protein = nan(4,13);
r2_p = nan(4,3); rmse_p = nan(4,3);

datasets = {'Feu_mean','Energy_mean','Protein_mean'};
impall_all = {impall_feu, impall_energy, impall_protein};
pdall_all = {}; xall_all = {};

for v = 1:3
    yielddata = eval(datasets{v});
    for p = 1:4
        disp(['Period ',num2str(p),' - ',datasets{v}])
        pp = Period_range(p,:) + 1; % Year indices
        
        % Build input data
        inputall = nan(3109,14,10);
        for i = 1:10
            yeard = pp(i);
            input = [precall(:,yeard), EDD_sum(:,yeard), GDD_sum(:,yeard), ...
                     notillall_int(:,yeard), irr_frac_int(:,yeard), nfer_frac_int(:,yeard), ...
                     tileall_frac3(:,yeard), ccall_frac3(:,yeard), ...
                     bdall, clayall, sandall, siltall, phall, yielddata(:,20+yeard)];
            inputall(:,:,i) = input;
        end
        input_m = nanmean(inputall,3);
        inputtable = array2table(input_m,'VariableNames', ...
            {'PRE','EDD','GDD','NT','IRR','NFER','TILE','CC','bd','clay','sand','silt','ph','YIELD'});
        
        [trainedModel,~] = RF_Model(inputtable,predictorNames);

        for j = 1:length(predictorNames)
            [pd,x] = partialDependence(trainedModel.RegressionEnsemble,predictorNames{j});
            pdall_all{v}(:,j,p) = pd'; xall_all{v}(:,j,p) = x';
        end

        res = trainedModel.predictFcn(inputtable(:,1:end-1));
        lm = fitlm(res,inputtable.YIELD);
        r2_p(p,v) = lm.Rsquared.Ordinary;
        rmse_p(p,v) = lm.RMSE;
        impall_all{v}(p,:) = predictorImportance(trainedModel.RegressionEnsemble);

        switch v
            case 1, Feu_p(:,p) = inputtable.YIELD;
            case 2, Energy_p(:,p) = inputtable.YIELD;
            case 3, Protein_p(:,p) = inputtable.YIELD;
        end
    end
end

contri_feu     = impall_all{1} ./ sum(impall_all{1},2);
contri_energy  = impall_all{2} ./ sum(impall_all{2},2);
contri_protein = impall_all{3} ./ sum(impall_all{3},2);


%% ===============================================================
%   STEP 3. Time Trends (1980–2020)
% ===============================================================

x = 1980:2020;
yields = {Feu_mean, Protein_mean, Energy_mean};
colors = {[88 177 103]/255,[66 133 244]/255,[234 97 35]/255};
titles = {'Food Equivalent Units','Protein','Energy'};
ylabs = {'FEU (kg/ha)','Protein yield (kg/ha)','Energy yield (MJ/ha)'};
contris = {contri_feu, contri_protein, contri_energy};

for v = 1:3
    y_m = nanmean(yields{v}(:,21:61));
    y_s = nanstd(yields{v}(:,21:61));
    figure('Color','w','Position',[563 278 784 497])
    hold on
    fill([x fliplr(x)], [y_m+y_s fliplr(y_m-y_s)], colors{v}, ...
        'FaceAlpha',0.3,'EdgeColor','none');
    plot(x, y_m, 'Color', colors{v}, 'LineWidth', 2);
    xlabel('Year'); ylabel(ylabs{v});
    title([titles{v},' Trend'],'FontSize',18);
    legend('± SD','Mean','Location','northwest','Box','off'); box on;
    % exportgraphics(gcf, sprintf('.\\Fig\\Fig_trend_%s.jpg',titles{v}), 'Resolution', 300)
   
    % --- Variable contribution subplot (stacked bar + pie)
    years = {'1980s','1990s','2000s','2010s'};
    stackData = contris{v}'*100;
    pieData = [sum(contris{v}(:,1:3),2),sum(contris{v}(:,4:8),2),sum(contris{v}(:,9:13),2)];
    barColors = [brewermap(3,'Oranges');brewermap(5,'Greens');brewermap(5,'Blues')];
    pieColors = [mean(brewermap(3,'Oranges'));mean(brewermap(5,'Greens'));mean(brewermap(5,'Blues'))];
    plotStackedBarWithTruePies(years, stackData, pieData, barColors, pieColors);
    % exportgraphics(gcf, sprintf('.\\Fig\\Fig_Contribution_%s.jpg',titles{v}), 'Resolution', 300)


end


%% ===============================================================
%   STEP 4. Spatial Patterns of Production by Decade (Enhanced)
% ===============================================================
%  This section visualizes spatial distributions of FEU, Energy, and Protein
%  for four decades (1980s–2010s). Each map includes:
%     - Decadal mean values per grid cell
%     - Colorbar with fixed range across decades
%     - An inset bar chart summarizing regional means
%     - Textbox showing statistical metrics (mean, R², RMSE)
%     - Variable category arrows indicating main drivers
% ===============================================================

regionShort = {'NE','SE','MW','GPN','GPS','SW','NW'};
softBlue = [198,219,239]/255;   % Soft blue for regional inset bars

% ---- Define value ranges across all decades for consistent color scaling ----
feu_range = [min(Feu_p(:)),    max(Feu_p(:))];
pro_range = [min(Protein_p(:)),max(Protein_p(:))];
ene_range = [min(Energy_p(:)), max(Energy_p(:))];

decades_names={'1980s','1990s','2000s','2010s'};
%% ===============================================================
%   FEU MAPS
% ===============================================================
for d = 1:4
    figure;
    color_value = Feu_p(:,d);

    % --- Base map ---
    ax = geoaxes;
    geoplot(USgeo.Shape, ColorData = color_value);
    ax.Basemap = 'none';
    ax.Grid = 'off';
    ax.Scalebar.Visible = 'off';
    ax.LongitudeAxis.Color = 'none';
    ax.LatitudeAxis.Color = 'none';
    ax.Position = [0.2, 0.25, 0.7, 0.55];

    % --- Colormap and colorbar ---
    mycolor = brewermap(100,"YlGn");
    colormap(mycolor);
    h = colorbar;
    h.Limits = feu_range;
    h.Label.String = 'FEU (kg/ha)';
    h.Label.FontSize = 14;
    h.FontSize = 14;
    h.Position = [0.81 0.3 0.02 0.22];
    box on;

    % --- Add arrows showing factor contributions (climate / management / soil) ---
    vars = struct;
    vars.top = {'PRE', 'EDD', 'GDD'};
    vars.top_contr = contri_feu(d,1:3)*100;
    vars.left = {'NT', 'IRR', 'NFER', 'DRA', 'CC'};
    vars.left_contr = contri_feu(d,4:8)*100;
    vars.bottom = {'BD', 'CLAY', 'SAND', 'SILT', 'PH'};
    vars.bottom_contr = contri_feu(d,9:13)*100;
    addVariableArrows(ax, vars);

    % --- Compute region-level statistics for inset bar chart ---
    data = array2table([UStable.RegionalID,color_value], 'VariableNames', {'id','value'});
    data_t = groupsummary(data,'id',@nanmean);

    % --- Model performance annotation ---
    R2 = r2_p(d,1);
    RMSE = rmse_p(d,1);
    meanvalue = nanmean(color_value);
    NRMSE = RMSE / meanvalue;
    infoText = sprintf('Mean = %2.2f\nR^2 = %.2f\nRMSE = %2.2f\n%%RMSE = %2.2f%%', ...
                       meanvalue, R2, RMSE, NRMSE*100);
    annotation('textbox', [0.68, 0.7, 0.2, 0.1], ...
        'String', infoText, 'FontSize', 11, 'FontWeight', 'normal', ...
        'EdgeColor', 'none', 'BackgroundColor', 'none', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    % --- Regional inset bar chart ---
    colors = mycolor(round(interp1(linspace(h.Limits(1),h.Limits(2),100),1:100,data_t.fun1_value)),:);
    insetPos = [0.23, 0.28, 0.2, 0.10];
    axInset = axes('Position', insetPos);
    b = bar(axInset, data_t.fun1_value, 'FaceColor', softBlue);
    b.FaceColor = 'flat'; b.CData = colors;
    axInset.XAxis.FontSize = 8;
    axInset.YAxis.FontSize = 8;
    xticks(1:numel(regionShort)); xticklabels(regionShort);
    set(gcf,'Position',[343 168 1059 697]);
    % exportgraphics(gcf, sprintf('.\\Fig\\Fig_FEU_MAP_%s.jpg',decades_names{d}), 'Resolution', 300)

end

%% ===============================================================
%   ENERGY MAPS
% ===============================================================
for d = 1:4
    figure;
    color_value = Energy_p(:,d) * 10 / 1000;

    ax = geoaxes;
    geoplot(USgeo.Shape, ColorData = color_value);
    ax.Basemap = 'none'; ax.Grid = 'off';
    ax.Scalebar.Visible = 'off';
    ax.LongitudeAxis.Color = 'none';
    ax.LatitudeAxis.Color = 'none';
    ax.Position = [0.2, 0.25, 0.7, 0.55];

    mycolor = brewermap(100,"YlGnBu");
    colormap(mycolor);
    h = colorbar;
    h.Limits = ene_range * 10 / 1000;
    h.Label.String = 'Energy (MJ/ha)';
    h.Label.FontSize = 14;
    h.Position = [0.81 0.3 0.02 0.22];
    box on;

    vars = struct;
    vars.top = {'PRE', 'EDD', 'GDD'};
    vars.top_contr = contri_energy(d,1:3)*100;
    vars.left = {'NT', 'IRR', 'NFER', 'DRA', 'CC'};
    vars.left_contr = contri_energy(d,4:8)*100;
    vars.bottom = {'BD', 'CLAY', 'SAND', 'SILT', 'PH'};
    vars.bottom_contr = contri_energy(d,9:13)*100;
    addVariableArrows(ax, vars);

    data = array2table([UStable.RegionalID,color_value],'VariableNames',{'id','value'});
    data_t = groupsummary(data,'id',@nanmean);

    R2 = r2_p(d,2);
    RMSE = rmse_p(d,2)*10/1000;
    meanvalue = nanmean(color_value);
    NRMSE = RMSE / meanvalue;
    infoText = sprintf('Mean = %2.2f\nR^2 = %.2f\nRMSE = %2.2f\n%%RMSE = %2.2f%%', ...
                       meanvalue,R2,RMSE,NRMSE*100);
    annotation('textbox',[0.68,0.7,0.2,0.1],'String',infoText,...
        'FontSize',11,'FontWeight','normal','EdgeColor','none','BackgroundColor','none');

    colors = mycolor(round(interp1(linspace(h.Limits(1),h.Limits(2),100),1:100,data_t.fun1_value)),:);
    insetPos = [0.23, 0.28, 0.2, 0.10];
    axInset = axes('Position', insetPos);
    b = bar(axInset, data_t.fun1_value, 'FaceColor', softBlue);
    b.FaceColor = 'flat'; b.CData = colors;
    axInset.XAxis.FontSize = 8; axInset.YAxis.FontSize = 8;
    xticks(1:numel(regionShort)); xticklabels(regionShort);
    set(gcf,'Position',[343 168 1059 697]);
     % exportgraphics(gcf, sprintf('.\\Fig\\Fig_Energy_MAP_%s.jpg',decades_names{d}), 'Resolution', 300)

end

%% ===============================================================
%   PROTEIN MAPS
% ===============================================================
for d = 1:4
    figure;
    color_value = Protein_p(:,d) * 10 / 1000;

    ax = geoaxes;
    geoplot(USgeo.Shape, ColorData = color_value);
    ax.Basemap = 'none'; ax.Grid = 'off';
    ax.Scalebar.Visible = 'off';
    ax.LongitudeAxis.Color = 'none';
    ax.LatitudeAxis.Color = 'none';
    ax.Position = [0.2, 0.25, 0.7, 0.55];

    mycolor = brewermap(100,"YlOrBr");
    colormap(mycolor);
    h = colorbar;
    h.Limits = pro_range * 10 / 1000;
    h.Label.String = 'Protein (kg/ha)';
    h.Label.FontSize = 14;
    h.Position = [0.81 0.3 0.02 0.22];
    box on;

    vars = struct;
    vars.top = {'PRE','EDD','GDD'};
    vars.top_contr = contri_protein(d,1:3)*100;
    vars.left = {'NT','IRR','NFER','DRA','CC'};
    vars.left_contr = contri_protein(d,4:8)*100;
    vars.bottom = {'BD','CLAY','SAND','SILT','PH'};
    vars.bottom_contr = contri_protein(d,9:13)*100;
    addVariableArrows(ax, vars);

    data = array2table([UStable.RegionalID,color_value],'VariableNames',{'id','value'});
    data_t = groupsummary(data,'id',@nanmean);

    R2 = r2_p(d,3);
    RMSE = rmse_p(d,3)*10/1000;
    meanvalue = nanmean(color_value);
    NRMSE = RMSE / meanvalue;
    infoText = sprintf('Mean = %2.2f\nR^2 = %.2f\nRMSE = %2.2f\n%%RMSE = %2.2f%%', ...
                       meanvalue,R2,RMSE,NRMSE*100);
    annotation('textbox',[0.68,0.7,0.2,0.1],'String',infoText,...
        'FontSize',11,'FontWeight','normal','EdgeColor','none','BackgroundColor','none');

    colors = mycolor(round(interp1(linspace(h.Limits(1),h.Limits(2),100),1:100,data_t.fun1_value)),:);
    insetPos = [0.23, 0.28, 0.2, 0.10];
    axInset = axes('Position', insetPos);
    b = bar(axInset, data_t.fun1_value, 'FaceColor', softBlue);
    b.FaceColor = 'flat'; b.CData = colors;
    axInset.XAxis.FontSize = 8; axInset.YAxis.FontSize = 8;
    xticks(1:numel(regionShort)); xticklabels(regionShort);
    set(gcf,'Position',[343 168 1059 697]);
     % exportgraphics(gcf, sprintf('.\\Fig\\Fig_Protein_MAP_%s.jpg',decades_names{d}), 'Resolution', 300)

end

%% ===============================================================
%   STEP 5. Partial Dependence Plots (PDP)
% ===============================================================
%  This section visualizes how each predictor variable (climate,
%  management, and soil) influences FEU, Energy, and Protein yields.
%  Each subplot shows a smoothed partial dependence curve for each
%  decade (1980s–2010s).
% ===============================================================

% --- Rename predictors with formatted labels ---
predictorNames2 = predictorNames;
predictorNames2{1}  = 'PRE (mm)';           % Precipitation
predictorNames2{2}  = 'EDD (°C days)';      % Extreme degree days
predictorNames2{3}  = 'GDD (°C days)';      % Growing degree days
predictorNames2{4}  = 'NT (%)';             % No-tillage fraction
predictorNames2{5}  = 'IRR (%)';            % Irrigation fraction
predictorNames2{6}  = 'NFER (kg N/ha)';     % Nitrogen fertilization
predictorNames2{7}  = 'DRA (%)';            % Drainage fraction
predictorNames2{8}  = 'CC (%)';             % Cover crop fraction
predictorNames2{9}  = 'BD (g/cm³)';         % Bulk density
predictorNames2{10} = 'CLAY (%)';
predictorNames2{11} = 'SAND (%)';
predictorNames2{12} = 'SILT (%)';
predictorNames2{13} = 'pH';

% ---------------------------------------------------------------
%   (1) FEU PDP Plot
% ---------------------------------------------------------------
figure('Color','w');
set(gcf,'Position',[311 207 1106 731]);

xall_feu=xall_all{1,1};
pdall_feu=pdall_all{1,1};
for i = 1:13
    subplot(4,4,i);
    hold on;

    % --- Loop through decades (1=1980s ... 4=2010s) ---
    for j = 1:4
        [xdata, ydata] = deal(xall_feu(:,i,j), pdall_feu(:,i,j));

        % Convert to relative change (%)
        ydata_c = (ydata - ydata(1)) * 100 ./ ydata(1);
        ydata_s = smooth(ydata_c, 50, 'sgolay', 3);

        % Convert fractions to percentage (for management variables)
        if i >= 4 && i <= 8 && i ~= 6
            xdata = xdata * 100;
        end

        % Plot PDP curve
        plot(xdata, ydata_s, '-', 'LineWidth', 1.5);
    end

    % --- Subplot formatting ---
    title(predictorNames2{i}, 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Change in FEU (%)', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end

% --- Unified legend placement ---
legend({'1980s','1990s','2000s','2010s'}, ...
       'FontSize', 12, 'Orientation', 'horizontal', ...
       'Position', [0.36 0.16 0.42 0.04]);

sgtitle('Partial Dependence of FEU on Climate, Management, and Soil Factors', ...
        'FontSize', 16, 'FontWeight', 'bold');
% exportgraphics(gcf, '.\Fig\Fig_FEU_PDP.jpg', 'Resolution', 300)

% ---------------------------------------------------------------
%   (2) ENERGY PDP Plot
% ---------------------------------------------------------------
figure('Color','w');
set(gcf,'Position',[311 207 1106 731]);

xall_energy=xall_all{1,2};
pdall_energy=pdall_all{1,2};
for i = 1:13
    subplot(4,4,i);
    hold on;
    for j = 1:4
        [xdata, ydata] = deal(xall_energy(:,i,j), pdall_energy(:,i,j));
        ydata_c = (ydata - ydata(1)) * 100 ./ ydata(1);
        ydata_s = smooth(ydata_c, 50, 'sgolay', 3);
        if i >= 4 && i <= 8 && i ~= 6
            xdata = xdata * 100;
        end
        plot(xdata, ydata_s, '-', 'LineWidth', 1.5);
    end
    title(predictorNames2{i}, 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Change in Energy (%)', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end

legend({'1980s','1990s','2000s','2010s'}, ...
       'FontSize', 12, 'Orientation', 'horizontal', ...
       'Position', [0.36 0.16 0.42 0.04]);

sgtitle('Partial Dependence of Energy Yield on Climate, Management, and Soil Factors', ...
        'FontSize', 16, 'FontWeight', 'bold');
% exportgraphics(gcf, '.\Fig\Fig_Energy_PDP.jpg', 'Resolution', 300)
% ---------------------------------------------------------------
%   (3) PROTEIN PDP Plot
% ---------------------------------------------------------------
figure('Color','w');
set(gcf,'Position',[311 207 1106 731]);
xall_protein=xall_all{1,3};
pdall_protein=pdall_all{1,3};
for i = 1:13
    subplot(4,4,i);
    hold on;
    for j = 1:4
        [xdata, ydata] = deal(xall_protein(:,i,j), pdall_protein(:,i,j));
        ydata_c = (ydata - ydata(1)) * 100 ./ ydata(1);
        ydata_s = smooth(ydata_c, 50, 'sgolay', 3);
        if i >= 4 && i <= 8 && i ~= 6
            xdata = xdata * 100;
        end
        plot(xdata, ydata_s, '-', 'LineWidth', 1.5);
    end
    title(predictorNames2{i}, 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Change in Protein (%)', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end

legend({'1980s','1990s','2000s','2010s'}, ...
       'FontSize', 12, 'Orientation', 'horizontal', ...
       'Position', [0.36 0.16 0.42 0.04]);

sgtitle('Partial Dependence of Protein Yield on Climate, Management, and Soil Factors', ...
        'FontSize', 16, 'FontWeight', 'bold');

% exportgraphics(gcf, '.\Fig\Fig_Protein_PDP.jpg', 'Resolution', 300)



