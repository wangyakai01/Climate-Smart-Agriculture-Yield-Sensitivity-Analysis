%% ==============================================================
%   Analysis of Crop Yield Sensitivity to Drought/Wet Years
%   Under Combined Management: Tillage + Cover Cropping + Drainage
% ==============================================================



%% ==============================================================
%   PART 1: Analysis for Dry Years (SPEI < 0)
% ==============================================================

% Initialize data containers
data1_1 = [];   % Only tillage
data1p = [];
data1_1_1 = []; % + Cover cropping
data1p_1 = [];
data1_2_1 = []; % + Cover cropping + Drainage (DRA)
data1p_2 = [];
crops={'Corn'	'Soybean'	'Wheat'};
group_inx=[repmat(1,1,11),repmat(2,1,11),repmat(3,1,11)];
% Loop over crop types (e.g., Corn and Soybean)
for c = 1:2
    k = 1;
    kres = []; pres = [];
    kres_1 = []; pres_1 = [];
    kres_2 = []; pres_2 = [];
    
    crop_type = crops(c);
    till_type_name = {'NT','RT','CT'};  % NT = No-till, RT = Reduced-till, CT = Conventional-till
    
    % Loop over years (2014–2024)
    for yeard = 2014:2024
        for j = 1:3
            try
                till_type = till_type_name{j};
                
                % Select subset: dry years, given crop type and tillage system
                index = Siteyield.SPEI_m < 0 & ...
                        Siteyield.year == yeard & ...
                        contains(Siteyield.Croptype, crop_type) & ...
                        contains(Siteyield.Tillage_type, till_type);
                Siteyield_s = Siteyield(index,:);
                
                % --- (1) Only Tillage: low CC + low DRA
                index1 = Siteyield_s.cover_crop_frac <= median(Siteyield_s.cover_crop_frac) & ...
                         Siteyield_s.TILE_frac <= median(Siteyield_s.TILE_frac);
                sres = fitlm(Siteyield_s.SPEI_m(index1), Siteyield_s.yield_adj_change(index1));
                kres(k,j) = sres.Coefficients.Estimate(2);
                pres(k,j) = sres.Coefficients.pValue(2);
                
                % --- (2) + Cover Cropping: high CC + low DRA
                index1 = Siteyield_s.cover_crop_frac >= median(Siteyield_s.cover_crop_frac) & ...
                         Siteyield_s.TILE_frac <= median(Siteyield_s.TILE_frac);
                sres = fitlm(Siteyield_s.SPEI_m(index1), Siteyield_s.yield_adj_change(index1));
                kres_1(k,j) = sres.Coefficients.Estimate(2);
                pres_1(k,j) = sres.Coefficients.pValue(2);
                
                % --- (3) + CC + DRA: high CC + high DRA
                index1 = Siteyield_s.cover_crop_frac > median(Siteyield_s.cover_crop_frac) & ...
                         Siteyield_s.TILE_frac > median(Siteyield_s.TILE_frac);
                sres = fitlm(Siteyield_s.SPEI_m(index1), Siteyield_s.yield_adj_change(index1));
                kres_2(k,j) = sres.Coefficients.Estimate(2);
                pres_2(k,j) = sres.Coefficients.pValue(2);
            catch
                continue;
            end
        end
        k = k + 1;
    end

    % Store results for each crop
    data1_1(:,c)   = kres(:);
    data1p(:,c)    = pres(:);
    data1_1_1(:,c) = kres_1(:);
    data1p_1(:,c)  = pres_1(:);
    data1_2_1(:,c) = kres_2(:);
    data1p_2(:,c)  = pres_2(:);
end

% Clean unrealistic values
data1_1(abs(data1_1) > 80 | data1_1 == 0) = nan;
data1_1_1(abs(data1_1_1) > 80 | data1_1_1 == 0) = nan;
data1_2_1(abs(data1_2_1) > 80 | data1_2_1 == 0) = nan;


%% ==============================================================
%   PART 2: Analysis for Wet Years (SPEI > 0)
% ==============================================================

data1_2 = [];
data1p = [];
data1_1_2 = [];
data1p_1 = [];
data1_2_2 = [];
data1p_2 = [];

for c = 1:2
    k = 1;
    kres = []; pres = [];
    kres_1 = []; pres_1 = [];
    kres_2 = []; pres_2 = [];
    
    crop_type = crops(c);
    till_type_name = {'NT','RT','CT'};
    
    for yeard = 2014:2024
        for j = 1:3
            try
                till_type = till_type_name{j};
                
                % Select subset: wet years, given crop type and tillage system
                index = Siteyield.SPEI_m > 0 & ...
                        Siteyield.year == yeard & ...
                        contains(Siteyield.Croptype, crop_type) & ...
                        contains(Siteyield.Tillage_type, till_type);
                Siteyield_s = Siteyield(index,:);
                
                % --- (1) Only Tillage
                index1 = Siteyield_s.cover_crop_frac <= median(Siteyield_s.cover_crop_frac) & ...
                         Siteyield_s.TILE_frac <= median(Siteyield_s.TILE_frac);
                sres = fitlm(Siteyield_s.SPEI_m(index1), Siteyield_s.yield_adj_change(index1));
                kres(k,j) = sres.Coefficients.Estimate(2);
                pres(k,j) = sres.Coefficients.pValue(2);
                
                % --- (2) + Cover Cropping
                index1 = Siteyield_s.cover_crop_frac >= median(Siteyield_s.cover_crop_frac) & ...
                         Siteyield_s.TILE_frac <= median(Siteyield_s.TILE_frac);
                sres = fitlm(Siteyield_s.SPEI_m(index1), Siteyield_s.yield_adj_change(index1));
                kres_1(k,j) = sres.Coefficients.Estimate(2);
                pres_1(k,j) = sres.Coefficients.pValue(2);
                
                % --- (3) + CC + DRA
                index1 = Siteyield_s.cover_crop_frac > median(Siteyield_s.cover_crop_frac) & ...
                         Siteyield_s.TILE_frac > median(Siteyield_s.TILE_frac);
                sres = fitlm(Siteyield_s.SPEI_m(index1), Siteyield_s.yield_adj_change(index1));
                kres_2(k,j) = sres.Coefficients.Estimate(2);
                pres_2(k,j) = sres.Coefficients.pValue(2);
            catch
                continue;
            end
        end
        k = k + 1;
    end
    
    data1_2(:,c)   = kres(:);
    data1p(:,c)    = pres(:);
    data1_1_2(:,c) = kres_1(:);
    data1p_1(:,c)  = pres_1(:);
    data1_2_2(:,c) = kres_2(:);
    data1p_2(:,c)  = pres_2(:);
end

% Clean unrealistic values
data1_2(abs(data1_2) > 80 | data1_2 == 0) = nan;
data1_1_2(abs(data1_1_2) > 80 | data1_1_2 == 0) = nan;
data1_2_2(abs(data1_2_2) > 80 | data1_2_2 == 0) = nan;


%% ==============================================================
%   PART 3: Visualization – Violin Plots of Yield Sensitivity
% ==============================================================

crop_title = {'Corn', 'Soybean'};
figure
set(gcf, 'Position', [236 224 1035 681])

% ---------------- Corn: Dry Years ----------------
subplot(2,2,1)
h = daviolinplot([data1_1(:,1), data1_1_1(:,1), data1_2_1(:,1)], ...
    'groups', group_inx, 'violin', 'full', 'scatter', 0, ...
    'scattersize', 30, 'boxalpha', 0.5, ...
    'xtlabels', {'Only Tillage','+CC','+CC+DRA'}, ...
    'legend', {'NT','RT','CT'});
hold on
plot([0,4],[0,0],'r--','LineWidth',2)
box on
h.lg.Visible = "off";
ylabel('Yield Change (%) per SPEI','FontSize',12,'FontWeight','bold')
xlabel('Management','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
title('Corn (Dry years)','FontSize',16,'FontWeight','bold');

% ---------------- Corn: Wet Years ----------------
subplot(2,2,2)
h = daviolinplot([data1_2(:,1), data1_1_2(:,1), data1_2_2(:,1)], ...
    'groups', group_inx, 'violin', 'full', 'scatter', 0, ...
    'scattersize', 30, 'boxalpha', 0.5, ...
    'xtlabels', {'Only Tillage','+CC','+CC+DRA'}, ...
    'legend', {'NT','RT','CT'});
hold on
plot([0,4],[0,0],'r--','LineWidth',2)
box on
h.lg.NumColumns = 3;
h.lg.Location = "north";
h.lg.Visible = "on";
h.lg.String = {'NT','RT','CT'};
ylabel('Yield Change (%) per SPEI','FontSize',12,'FontWeight','bold')
xlabel('Management','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
title('Corn (Wet years)','FontSize',16,'FontWeight','bold');

% ---------------- Soybean: Dry Years ----------------
subplot(2,2,3)
h = daviolinplot([data1_1(:,2), data1_1_1(:,2), data1_2_1(:,2)], ...
    'groups', group_inx, 'violin', 'full', 'scatter', 0, ...
    'scattersize', 30, 'boxalpha', 0.5, ...
    'xtlabels', {'Only Tillage','+CC','+CC+DRA'}, ...
    'legend', {'NT','RT','CT'});
hold on
plot([0,4],[0,0],'r--','LineWidth',2)
box on
h.lg.Visible = "off";
ylabel('Yield Change (%) per SPEI','FontSize',12,'FontWeight','bold')
xlabel('Management','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
title('Soybean (Dry years)','FontSize',16,'FontWeight','bold');

% ---------------- Soybean: Wet Years ----------------
subplot(2,2,4)
h = daviolinplot([data1_2(:,2), data1_1_2(:,2), data1_2_2(:,2)], ...
    'groups', group_inx, 'violin', 'full', 'scatter', 0, ...
    'scattersize', 30, 'boxalpha', 0.5, ...
    'xtlabels', {'Only Tillage','+CC','+CC+DRA'}, ...
    'legend', {'NT','RT','CT'});
hold on
plot([0,4],[0,0],'r--','LineWidth',2)
box on
h.lg.Visible = "off";
ylabel('Yield Change (%) per SPEI','FontSize',12,'FontWeight','bold')
xlabel('Management','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold')
title('Soybean (Wet years)','FontSize',16,'FontWeight','bold');

% ---------------- Formatting ----------------
for i = 1:4
    subplot(2,2,i)
    ylim([-150,150])
    text(0.6,170,sprintf('%s', 'a'+i-1), 'FontSize',16,'FontWeight','bold')
end

% Save figure
 % exportgraphics(gcf, '.\Fig\Fig_site_Result_dry_wetyears_update.jpg', 'Resolution', 300)
