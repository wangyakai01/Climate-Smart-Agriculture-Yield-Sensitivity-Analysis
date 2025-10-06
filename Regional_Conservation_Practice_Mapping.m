%% ============================================================
%   Regional Conservation Practice Mapping (1980s–2010s)
%   Practices: No-Till, Tile Drainage, and Cover Cropping
% ============================================================
%  Notes:
%   - USgeo: shapefile containing U.S. county boundaries
%   - *_int / *_frac3: fractional area of practice (0–1) per county per year
%   - Decadal ranges:
%         2:11   → 1980s
%         12:21  → 1990s
%         22:31  → 2000s
%         32:41  → 2010s
% ============================================================
clear; clc;

%% ===============================================================
% Define practices, plotting configurations, and scaling
practices = {
    'No-Till',          notillall_int,  [0 80],  "Reds",   "-RdBu",  'nt';
    'Tile Drainage',    tileall_frac3,  [0 60],  "Oranges","-PuOr",  'tile';
    'Cover Cropping',   ccall_frac3,    [0 20],  "Greens", "PRGn",   'cc'
};

% Define decades
decades = {
    '1980s',  2:11;
    '1990s', 12:21;
    '2000s', 22:31;
    '2010s', 32:41
};

% Loop through each practice
for k = 1:size(practices,1)
    practice_name = practices{k,1};
    data_matrix   = practices{k,2};
    clim_base     = practices{k,3};
    cmap_base     = practices{k,4};
    cmap_change   = practices{k,5};
    short_name    = practices{k,6};

    fprintf('\n--- Generating maps for %s ---\n', practice_name);

    %% ========== Decadal Mean Maps ==========
    for d = 1:size(decades,1)
        dec_name = decades{d,1};
        year_range = decades{d,2};

        figure('Color','w');
        set(gcf,'Position',[408   426   899   464]);

        % Compute mean fraction (%)
        color_value = nanmean(data_matrix(:,year_range), 2) * 100;

        % ---- Base map ----
        ax = geoaxes;
        geoplot(USgeo.Shape, ColorData=color_value);
        ax.Basemap = 'none'; ax.Grid = 'off';
        ax.Scalebar.Visible = 'off';
        ax.LongitudeAxis.Color = 'none'; ax.LatitudeAxis.Color = 'none';

        % ---- Colormap and colorbar ----
        colormap(brewermap([], cmap_base));
        h = colorbar;
        h.Position = [0.8 0.15 0.0212 0.4];
        h.Label.String = sprintf('%s area (%%, %s)', practice_name, dec_name);
        h.Label.FontSize = 12;
        box off; clim(clim_base);
        title(sprintf('%s — %s', practice_name, dec_name), 'FontSize', 14, 'FontWeight', 'bold');

        % ---- Optional export ----
        %exportgraphics(gcf, sprintf('.\\Fig\\Fig_%s_%s.jpg', short_name, dec_name), 'Resolution', 300);
    end

    %% ========== Change Map (2010s − 1980s) ==========
    figure('Color','w');
    set(gcf,'Position',[680 569 573 309]);

    color_value = (nanmean(data_matrix(:,32:41),2) - nanmean(data_matrix(:,2:11),2)) * 100;

    ax = geoaxes;
    geoplot(USgeo.Shape, ColorData=color_value);
    ax.Basemap = 'none'; ax.Grid = 'off';
    ax.Scalebar.Visible = 'off';
    ax.LongitudeAxis.Color = 'none'; ax.LatitudeAxis.Color = 'none';

    colormap(brewermap([], cmap_change));
    h = colorbar;
    h.Position = [0.8 0.15 0.0212 0.4];
    h.Label.String = sprintf('%s area change (%%, 2010s−1980s)', practice_name);
    h.Label.FontSize = 12;
    box off; clim([-clim_base(2), clim_base(2)]);
    title(sprintf('%s — Change (2010s − 1980s)', practice_name), 'FontSize', 14, 'FontWeight', 'bold');

    %exportgraphics(gcf, sprintf('.\\Fig\\Fig_%s_change.jpg', short_name), 'Resolution', 300);
end

