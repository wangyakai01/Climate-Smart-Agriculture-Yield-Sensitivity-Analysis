%% ===============================================================
%     RANDOM FOREST ANALYSIS OF CROP YIELD RESPONSE TO MANAGEMENT
% ===============================================================
% This script:
% 1. Trains Random Forest (RF) models for each crop and time period.
% 2. Evaluates management effects (No-till, Cover Crop, Drainage, Irrigation).
% 3. Computes drought sensitivity (ADSI) across regions and crops.
% ===============================================================


%% ===============================================================
%   PART 1: Train Random Forest for Each Crop Type
% ===============================================================
rfresall_c = struct;
for i = 1:7
    m = keyperiod(i).keymonth;             % Key month for this crop
    yielddata = eval(cropyieldnames_f{i}); % Crop yield dataset

    % Preallocate storage for each year
    impall_c1 = nan(41, 13);
    validationRMSE1 = nan(41, 1);
    r21 = nan(41, 1);
    rmse1 = nan(41, 1);
    res = nan(3109, 41);          % Baseline yield prediction
    res1 = nan(3109, 41, 6);      % Predictions under management removal (zero)
    res2 = nan(3109, 41, 6);      % Predictions under management addition (one)

    % Climate and management data
    precall = metindexall(i).prec;
    GDD_sum = metindexall(i).gdd;
    EDD_sum = metindexall(i).edd;

    for yeard = 1:41
        disp([i, yeard])

        % Climate indicators
        spei = speidata(m).spei(:, yeard);

        % ==== Combine predictors ====
        input = [precall(:, yeard), EDD_sum(:, yeard), GDD_sum(:, yeard), spei, ...
                 notillall_int(:, yeard), irr_frac_int(:, yeard), nfer_frac_int(:, yeard), ...
                 tileall_frac3(:, yeard), ccall_frac3(:, yeard),noydata_int(:,yeard),nhxdata_int(:,yeard), ...
                 bdall, clayall, sandall, siltall, phall, yielddata(:, 20 + yeard)];

        inputtable = array2table(input, ...
            'VariableNames', {'PRE','EDD','GDD','SPEI','NT','IRR','NFER','TILE','CC','NOY','NHY','bd','clay','sand','silt','ph','YIELD'});

        predictorNames = {'PRE','EDD','GDD','NT','IRR','NFER','TILE','CC','bd','clay','sand','silt','ph'};

        % ==== Train Random Forest Model ====
        [trainedModel, validationRMSE1(yeard)] = RF_Model(inputtable, predictorNames);

        % ==== Predict baseline yield ====
        inputdata = inputtable(:, 1:end-1);
        res(:, yeard) = trainedModel.predictFcn(inputdata);

        % ==== Management effect simulations (set 0 = remove, 1 = add) ====
        inputdata_o = inputdata;

        % --- Remove individual or combined managements ---
        configs_zero = {
            {'IRR'}, {'NT'}, {'TILE'}, {'CC'}, ...
            {'NT','CC','TILE'}, {'NT','CC'}
        };
        for s = 1:6
            inputdata = inputdata_o;
            for f = 1:length(configs_zero{s})
                
                inputdata.(configs_zero{s}{f}) = zeros(3109,1);
            end
            res1(:, yeard, s) = trainedModel.predictFcn(inputdata);
        end

        % --- Add managements (set to 1) ---
        configs_one = {
            {'IRR'}, {'NT'}, {'TILE'}, {'CC'}, ...
            {'NT','CC','TILE'}, {'NT','CC'}
        };
        for s = 1:6
            inputdata = inputdata_o;
            for f = 1:length(configs_one{s})
                inputdata.(configs_one{s}{f}) = ones(3109,1);
            end
            res2(:, yeard, s) = trainedModel.predictFcn(inputdata);
        end

        % ==== Evaluate performance ====
        lm = fitlm(res(:, yeard), inputtable.YIELD);
        r21(yeard) = lm.Rsquared.Ordinary;
        rmse1(yeard) = lm.RMSE;

        % ==== Predictor importance ====
        imp = predictorImportance(trainedModel.RegressionEnsemble);
        impall_c1(yeard, :) = imp;
    end

    % Store model results for this crop
    rfresall_c(i).r2 = r21;
    rfresall_c(i).rmse = rmse1;
    rfresall_c(i).res = res;
    rfresall_c(i).res1 = res1;
    rfresall_c(i).res2 = res2;
    rfresall_c(i).impall_c1 = impall_c1;
    rfresall_c(i).validationRMSE1 = validationRMSE1;
end

    
%% ===============================================================
%   PART 2: Compute Yield–SPEI Sensitivity (with vs. without CSA)
% ===============================================================
for type = 1:7
    disp(type)
    sen_spei_m = crop_yield_spei(type).sen_spei_m;
    res = rfresall_c(type).res;
   
    res2 = rfresall_c(type).res2;

    corn_yield_kall = struct;
    spei_s = crop_yield_spei(type).speivalues;

    for dd = 1:6  % 6 management configurations
        noncsa = res;
        csa = res2(:,:,dd);

        
        corn_yield_b1 = nan(3109,2);
        corn_yield_k1 = nan(3109,2);

        warning off
        for k = 1:3109
            % disp([dd, k])
            data1 = noncsa(k,:);
            data2 = csa(k,:);
            data3 = spei_s(k,1:41);

            if all(isnan(data3)), continue, end

            % Remove long-term trend (detrend using Savitzky-Golay)
            data1_climate = (data1 - smooth(data1,41,'sgolay',3)') ./ smooth(data1,41,'sgolay',3)';
            data2_climate = (data2 - smooth(data2,41,'sgolay',3)') ./ smooth(data2,41,'sgolay',3)';

            p = fitlm(data3', data1_climate');
            corn_yield_b1(k,1) = p.Coefficients.Estimate(1);
            corn_yield_k1(k,1) = p.Coefficients.Estimate(2);

            p = fitlm(data3', data2_climate');
            corn_yield_b1(k,2) = p.Coefficients.Estimate(1);
            corn_yield_k1(k,2) = p.Coefficients.Estimate(2);
        end

        corn_yield_kall(dd).data = corn_yield_k1;
    end
    crop_yield_spei(type).yield_kall = corn_yield_kall;
end


%% ===============================================================
%   PART 3: Compute Agricultural Drought Sensitivity Index (ADSI)
% ===============================================================
crop_yield_k2 = [];
crop_yield_k1 = [];

for type = 1:7
    corn_yield_kall = crop_yield_spei(type).yield_kall;
    kk = 1;
    for dd = [2,6,5] % Management configs of interest
        corn_yield_k1 = corn_yield_kall(dd).data;
        crop_yield_k2(:,type,kk) = corn_yield_k1(:,2); % CSA scenario
        crop_yield_k1(:,type,kk) = corn_yield_k1(:,1); % Non-CSA
        kk = kk + 1;
    end
end


% ==== Area-weighted ADSI aggregation ====
areanames = {'cornareadata1','soybeanareadata1','wheat_winterareadata1', ...
             'wheat_springareadata1','SORGHUMareadata1','BARLEYareadata1','RICEareadata1'};
frac1 = [];
for i = 1:7
    cropfrac = eval(areanames{i});
    frac1(:,i) = nanmean(cropfrac(:,21:61), 2);
end

% Area-weighted means across crops
 data1=crop_yield_k1(:,:,1);
data11=nansum(data1.*frac1,2)./nansum(frac1,2);
data11(data11==0)=nan;
data2=crop_yield_k2(:,:,1);
data21=nansum(data2.*frac1,2)./nansum(frac1,2);
data21(data21==0)=nan;
data2=crop_yield_k2(:,:,2);
data22=nansum(data2.*frac1,2)./nansum(frac1,2);
data22(data22==0)=nan;
data2=crop_yield_k2(:,:,3);
data23=nansum(data2.*frac1,2)./nansum(frac1,2);
data23(data23==0)=nan;


% ==== Compute management effects on ADSI ====
eff(:,1) = data21 - data11;
eff(:,2) = data22 - data11;
eff(:,3) = data23 - data11;

data = array2table([UStable.RegionalID, eff*100], 'VariableNames', {'id','value1','value2','value3'});
data_t = groupsummary(data, 'id', @nanmean);


%% ===============================================================
%   PART 4: Spatial Mapping of ADSI and Management Effects
% ===============================================================

figure
proj = projcrs(102008, Authority="ESRI");
ax = newmap(proj);

color_value = data11 * 100;
geoplot(ax, USgeo.Shape, ColorData=color_value, EdgeColor=[0.7 0.7 0.7]);
hold on
for r = 1:7
    geoplot(ax, USregion(r).Y, USregion(r).X, '-', "LineWidth", 2, 'Color', [0.5,0.5,0.5])
end

ax.Scalebar.Visible = 'off';
ax.Visible = "off";
colormap(brewermap([], 'YlOrRd'));
set(gcf, 'Position', [1 49 1920 955])

cl = colorbar;
cl.Location = 'south';
cl.Position = [0.24 0.06 0.2 0.0223];
cl.Label.String = {'Agricultural Drought Sensitivity Index (ADSI)'};
cl.Label.Position = [7.4 2 0];
cl.FontSize = 12;

% ==== Regional inset bar charts ====
regionShort = {'NE','SE','MW','GPN','GPS','SW','NW','US'};
insetposall = [
    0.76, 0.53, 0.05, 0.1;
    0.72, 0.28, 0.05, 0.1;
    0.65, 0.72, 0.05, 0.1;
    0.422, 0.12, 0.05, 0.1;
    0.52, 0.12, 0.05, 0.1;
    0.355, 0.18, 0.05, 0.1;
    0.28, 0.25, 0.05, 0.1
];

for ii = 1:7
    if isnan(data_t{ii,3}), continue, end
    axInset = axes('Position', insetposall(ii,:));
    b = bar(axInset, data_t{ii,3:end}, 'FaceAlpha', 0.5);
    b.FaceColor = 'flat';
    b.CData = brewermap(3, 'paired');
    ylabel('Change in ADSI')
    set(gca, 'Color', 'none', 'Box', 'off', 'XColor', 'none')
    title(regionShort{ii}, "FontSize", 12)
end

% ==== National summary inset ====
eff_m = nanmean(eff) * 100;
mycolors = brewermap(3, 'paired');
labels = {'NT', 'NT+CC', 'NT+CC+DRA'};
axInset = axes('Position', [0.73, 0.12, 0.05, 0.1]);
hold on
for f = 1:3
    bar(axInset, f, eff_m(f), 'FaceColor', mycolors(f,:), 'FaceAlpha', 0.5, 'DisplayName', labels{f});
end
ylabel('Change in ADSI')
set(gca, 'Color', 'none', 'Box', 'off', 'XColor', 'none')
title(regionShort{8}, "FontSize", 12)

% Legend
le = legend;
le.Box = 'off';
le.Orientation = "vertical";
le.Position = [0.22, 0.15, 0.1344, 0.035];
le.FontSize = 12;
%%
exportgraphics(gcf,sprintf('.\\Fig\\Fig_csa_%s.jpg','ADSI'),'Resolution',300)

