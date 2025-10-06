%% ========================== Figure 2 ==========================


regionShort = {'NE','SE','MW','GPN','GPS','SW','NW'}; 
%% ========================== Data Preparation ==========================
% List of datasets: crop yield and nutrient data
crop_names = {
    'cornyielddata','soybeanyielddata','wheat_winteryielddata', ...
    'wheat_springyielddata','SORGHUMyielddata','BARLEYyielddata', ...
    'RICEyielddata','Feu_mean','Energy_mean','Protein_mean'
};

% Initialize structure to store correlation and regression results
crop_yield_corr = struct;

for ii = 1:length(crop_names)
    
    cropyielddata = eval(crop_names{ii});

    % Preallocate matrices for correlation, slope (k), and intercept (b)
    num_counties = 3109;
    num_months = 12;
    nan_matrix = nan(num_counties, num_months);

    % Create 7 time segments for analysis
    for seg = 1:7
        eval(sprintf('corn_yield_corr%d = nan_matrix;', seg));
        eval(sprintf('corn_yield_k%d = nan_matrix;', seg));
        eval(sprintf('corn_yield_b%d = nan_matrix;', seg));
    end

    % Loop through all counties
    for i = 1:height(cropyielddata)
        warning off
disp([ii,i])
        % Detrend yield data using a Savitzky-Golay filter (63-point window)
        yield_climate = (cropyielddata(i,:) - smooth(cropyielddata(i,:),63,'sgolay',3)') ...
            ./ smooth(cropyielddata(i,:),63,'sgolay',3)';
        yield_climate1 = yield_climate(21:end)';

        % ===== Correlation windows =====
        % Each block analyzes a specific time window and computes
        % correlation + linear regression between yield anomalies and SPEI
        windows = {
            [1,40],  'corr1', [1,40]; ...
            [1,20],  'corr2', [1,20]; ...
            [21,40], 'corr3', [21,40]; ...
            [1,10],  'corr4', [1,10]; ...
            [11,20], 'corr5', [11,20]; ...
            [21,30], 'corr6', [21,30]; ...
            [31,40], 'corr7', [31,40]
        };

        % Loop over different time windows
        for w = 1:size(windows,1)
            idx = windows{w,3};
            if sum(~isnan(yield_climate1(idx(1):idx(2)))) > 5
                for j = 1:num_months
                    spei_s = speidata(j).spei(i,:);
                    [xdata, ydata] = prepareCurveData(spei_s(idx(1):idx(2))', yield_climate1(idx(1):idx(2)));
                    corr_val = corr(xdata, ydata);
                    p = fitlm(xdata, ydata);
                    eval(sprintf('corn_yield_corr%d(i,j) = corr_val;', w));
                    eval(sprintf('corn_yield_b%d(i,j) = p.Coefficients.Estimate(1);', w));
                    eval(sprintf('corn_yield_k%d(i,j) = p.Coefficients.Estimate(2);', w));
                end
            end
        end
    end

    % ===== Store results into structure =====
    for seg = 1:7
        crop_yield_corr(ii).(sprintf('corr%d', seg)) = eval(sprintf('corn_yield_corr%d', seg));
        crop_yield_corr(ii).(sprintf('k%d', seg)) = eval(sprintf('corn_yield_k%d', seg));
        crop_yield_corr(ii).(sprintf('b%d', seg)) = eval(sprintf('corn_yield_b%d', seg));
    end
end

%% ========================== Figure: FEU 1980s ==========================
figure
set(gcf, 'Position', [680 569 573 309])
color_value = crop_yield_corr(8).k4(:,8);  % FEU sensitivity (1980s)
ax = geoaxes;
geoplot(USgeo.Shape, ColorData=color_value);
ax.Basemap = 'none';
ax.Grid = 'off';
ax.Scalebar.Visible = 'off';
ax.LongitudeAxis.Color = 'none';
ax.LatitudeAxis.Color = 'none';
colormap(brewermap([],"YlOrRd"))
h = colorbar;
h.Position = [0.8 0.15 0.0212 0.4];
box off
clim([-0.15, 0.3])
%exportgraphics(gcf, '.\Fig\Fig_FEU_1980s.jpg', 'Resolution', 300)

%% ========================== Figure: FEU 2010s ==========================
figure
set(gcf, 'Position', [680 569 573 309])
color_value = crop_yield_corr(8).k7(:,8);  % FEU sensitivity (2010s)
ax = geoaxes;
geoplot(USgeo.Shape, ColorData=color_value);
ax.Basemap = 'none';
ax.Grid = 'off';
ax.Scalebar.Visible = 'off';
ax.LongitudeAxis.Color = 'none';
ax.LatitudeAxis.Color = 'none';
colormap(brewermap([],"YlOrRd"))
h = colorbar;
h.Position = [0.8 0.15 0.0212 0.4];
box off
clim([-0.15, 0.3])
%exportgraphics(gcf, '.\Fig\Fig_FEU_2010s.jpg', 'Resolution', 300)

%% ========================== Figure: FEU Change ==========================
figure
set(gcf, 'Position', [680 569 573 309])
color_value = crop_yield_corr(8).k7(:,8) - crop_yield_corr(8).k4(:,8); % Change between decades
ax = geoaxes;
geoplot(USgeo.Shape, ColorData=color_value);
ax.Basemap = 'none';
ax.Grid = 'off';
ax.Scalebar.Visible = 'off';
ax.LongitudeAxis.Color = 'none';
ax.LatitudeAxis.Color = 'none';
colormap(brewermap([],"-RdYlBu"))
h = colorbar;
h.Position = [0.8 0.15 0.0212 0.4];
box off
clim([-0.3, 0.3])
%exportgraphics(gcf, '.\Fig\Fig_FEU_s.jpg', 'Resolution', 300)

%% ========================== Regional Statistics ==========================
color_value = crop_yield_corr(8).k7(:,8) - crop_yield_corr(8).k4(:,8);
data = array2table([UStable.RegionalID, color_value], 'VariableNames', {'id','value'});

% Compute mean and standard deviation by region
data_s = groupsummary(data, 'id', @nanstd);
data_t = groupsummary(data, 'id', @nanmean);
data_t.std = data_s.fun1_value;
data_t.regionShort = regionShort';
data_t = sortrows(data_t, "fun1_value", "ascend");

% Map colors based on value
mycolor = brewermap([],"-RdYlBu");
colors = mycolor(round(interp1(linspace(-0.3,0.3,256),1:256,data_t.fun1_value)),:);

% ===== Plot regional mean change =====
figure
set(gcf, 'Position', [680 569 573 309])
b = bar(1:7, data_t.fun1_value, 'FaceColor', 'flat');
hold on

% Add error bars
errors = data_t.std;
errorbar(1:7, data_t.fun1_value, errors, errors, 'LineWidth', 1.5, ...
    'Color', 'k', 'CapSize', 12, 'LineStyle', 'none');

b.CData = colors;
ax = gca;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
xticks(1:7);
xticklabels(data_t.regionShort);
ylabel('Change in yield sensitivity to SPEI')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
%exportgraphics(gcf, '.\Fig\Fig_FEU_change.jpg', 'Resolution', 300)

