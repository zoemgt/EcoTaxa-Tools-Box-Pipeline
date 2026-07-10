%% EcoTaxa Tool Box — main pipeline (MATLAB 2025b compatible)
%
% Processes EcoTaxa exports from quantitative imaging instruments
% (UVP, IFCB, FlowCam, PlanktoScope, ZooScan) and produces ecological
% figures: maps, abundance/biovolume bar charts, NBSS spectra,
% and multivariate analyses.
%
% Reference pipeline: github.com/zoemgt/ecotaxatools

clear all
close all
warning('off');

% Add lib/ to path
currentFolder = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(currentFolder, 'lib')));
[directoryoftoolbox, ~, ~] = fileparts(mfilename('fullpath'));

%% ── 1. Data loading ─────────────────────────────────────────────────────

list = {'EcoTaxa export (tsv files) — process now', 'Matlab base.mat — already processed'};
[datatoload, tf] = listdlg('PromptString', 'What data should be loaded?', ...
    'ListString', list, 'SelectionMode', 'single', 'ListSize', [400, 200]);
if ~tf, return; end

if datatoload == 1
    ecotaxa_tools;  % runs TSV processing -> produces base in workspace

    list2 = {'Stop after processing', 'Continue to ecological analysis'};
    [indx, tf] = listdlg('PromptString', 'Next step?', ...
        'ListString', list2, 'SelectionMode', 'single', 'ListSize', [400, 200]);
    if ~tf || indx == 1, return; end

    % Reload the saved base for the analysis section
    clear all
    f = msgbox('Select the saved *.mat base file');
    uiwait(f);
    [file, path] = uigetfile('*.mat');
    load(fullfile(path, file));

else
    f = msgbox('Select the saved *.mat base file');
    uiwait(f);
    [file, path] = uigetfile('*.mat');
    load(fullfile(path, file));
end

m = numel(base);

%% ── 2. Sort samples by date ─────────────────────────────────────────────

[~, sort_idx] = sort([base.date]);
base = base(sort_idx);

%% ── 3. Geographic coordinates and sample labels ──────────────────────────

lat = [base.latitude];
lon = [base.longitude];

sampleid = {};
for i = 1:m
    sampleid = [sampleid; base(i).sampleid];
end

%% ── 4. Export option ─────────────────────────────────────────────────────

list2 = {'Yes', 'No'};
[saveChoice, tf] = listdlg('PromptString', 'Export all graphs as PDF?', ...
    'ListString', list2, 'SelectionMode', 'single', 'ListSize', [400, 200]);
if ~tf, saveChoice = 2; end
doSave = (saveChoice == 1);

%% ── 5. Map: station locations ────────────────────────────────────────────

scrsz     = get(0, 'ScreenSize');
landcolor = [1, 0.949, 0.867];

latmax = max(lat);  latmin = min(lat);
lonmax = max(lon);  lonmin = min(lon);
lat_ext = max(0.5, 0.30 * (latmax - latmin));
lon_ext = max(0.5, 0.30 * (lonmax - lonmin));
latmax = latmax + lat_ext;  latmin = latmin - lat_ext;
lonmax = lonmax + lon_ext;  lonmin = lonmin - lon_ext;

figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
m_proj('mercator', 'lon', [lonmin lonmax], 'lat', [latmin latmax]);
m_coast('patch', landcolor, 'edgecolor', [0.6 0.6 0.6]);
m_grid('fontsize', 8, 'linestyle', 'none', 'meridianlabel', 'off', ...
       'parallellabel', 'off', 'xtick', [], 'ytick', []);
hold on
m_plot(lon, lat, 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.5 0.5 0.5]);
m_text(lon, lat, sampleid);
if doSave
    exportgraphics(gcf, 'Map of stations.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(2); close gcf

%% ── 6. Build abundance / biovolume matrices ──────────────────────────────

Ab = [];
Bv = [];

for i = 1:m
    Ab = [Ab; base(i).regroupped.Ab];
    Bv = [Bv; base(i).regroupped.Bv];
    plankton_groups = base(i).regroupped.Zoo_groups;
    nStart = size(base(1).regroupped.Idlistshort, 1) + 1;
    taxo_annotation_shortname = [base(i).regroupped.Idlistshort; ...
                                  base(i).regroupped.Zoo_groups(nStart:end)];
end

% livingplace is a 2-element range vector [start end]: use start for "all living" column
lp = base(1).regroupped.livingplace;
Ab_tot = Ab(:, lp(1));
Bv_tot = Bv(:, lp(1));

%% ── 7. Extract and filter analysis subsets ───────────────────────────────

% Taxonomic level
op = base(1).regroupped.originalplace;
data_ab_taxo          = Ab(:, op(1):op(2));
data_bv_taxo          = Bv(:, op(1):op(2));
plankton_groups_data_taxo = taxo_annotation_shortname(op(1):op(2));

[data_ab_taxo, data_bv_taxo, plankton_groups_data_taxo, excl_idx_taxo] = ...
    filter_groups(data_ab_taxo, data_bv_taxo, plankton_groups_data_taxo, 'taxonomic annotation');

% PFT group level
pp = base(1).regroupped.planktongroupplace;
data_ab_group          = Ab(:, pp(1):pp(2));
data_bv_group          = Bv(:, pp(1):pp(2));
plankton_groups_data_group = taxo_annotation_shortname(pp(1):pp(2));

% Detritus first by default
idx_det = find(strcmp(plankton_groups_data_group, 'detritus'));
if ~isempty(idx_det)
    col_order = [idx_det, setdiff(1:numel(plankton_groups_data_group), idx_det)];
    plankton_groups_data_group = plankton_groups_data_group(col_order);
    data_ab_group = data_ab_group(:, col_order);
    data_bv_group = data_bv_group(:, col_order);
else
    col_order = 1:numel(plankton_groups_data_group);
end

[data_ab_group, data_bv_group, plankton_groups_data_group, excl_idx_group] = ...
    filter_groups(data_ab_group, data_bv_group, plankton_groups_data_group, 'plankton groups');

% Trophic level
tp = base(1).regroupped.trophicplace;
data_ab_troph          = Ab(:, tp(1):tp(2));
data_bv_troph          = Bv(:, tp(1):tp(2));
plankton_groups_data_troph = taxo_annotation_shortname(tp(1):tp(2));

[data_ab_troph, data_bv_troph, plankton_groups_data_troph, excl_idx_troph] = ...
    filter_groups(data_ab_troph, data_bv_troph, plankton_groups_data_troph, 'trophic groups');

% Trophic readable labels
trophic_labels = { '-1','not feeding'; '1','phototrophs'; '1.5','mixotrophs'; ...
                    '2','grazers'; '2.5','omnivorous'; '3','predators'; '3.5','unknown' };
plankton_groups_data_troph_name = plankton_groups_data_troph;
for i = 1:numel(plankton_groups_data_troph)
    row = find(strcmp(strtrim(plankton_groups_data_troph{i}), trophic_labels(:,1)), 1);
    if ~isempty(row)
        plankton_groups_data_troph_name{i} = trophic_labels{row, 2};
    end
end

%% ── 8. Per-station community summary ─────────────────────────────────────

ESD = base(1).tot.ESDvector;

for i = 1:m

    figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
    p_group_colors;

    sid = strrep(char(base(i).sampleid), '_', ' ');

    % Shannon index (on taxonomic groups)
    ab_row = data_ab_taxo(i,:);
    ab_tot_i = sum(ab_row, 'omitnan');
    if ab_tot_i > 0
        relAb_taxo = ab_row ./ ab_tot_i;
        H = -sum(relAb_taxo .* log(relAb_taxo), 'omitnan');
    else
        H = NaN;
    end

    % Relative abundance per PFT group (for pie)
    ab_group_row = data_ab_group(i,:);
    ab_group_tot = sum(ab_group_row, 'omitnan');
    if ab_group_tot > 0
        relAb_group = 100 * ab_group_row ./ ab_group_tot;
    else
        relAb_group = zeros(size(ab_group_row));
    end

    % Panel 1 — station info
    subplot(3,4,1); axis off
    text(0, 0.8, sid, 'FontSize', 11, 'FontWeight', 'bold');
    text(0, 0.5, ['Total abundance (ind m^{-3}): '  num2str(Ab_tot(i), '%.2f')], 'FontSize', 10);
    text(0, 0.3, ['Total biovolume (mm^3 m^{-3}): ' num2str(Bv_tot(i), '%.4f')], 'FontSize', 10);
    text(0, 0.1, ['Shannon index: ' num2str(H, '%.2f')],                          'FontSize', 10);

    % Panel 2 — NBSS
    ax2 = subplot(3,4,2);
    pos2 = get(ax2, 'Position');
    set(ax2, 'Position', [pos2(1)+0.03, pos2(2), pos2(3), pos2(4)]);
    nbss_data = base(i).regroupped.Ybv_Ellipsoid_BV_spectra(:, base(i).regroupped.livingplace(1));
    valid_nbss = nbss_data > 0;
    if any(valid_nbss)
        plot(ESD(valid_nbss), nbss_data(valid_nbss), ...
             'o', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
        set(gca, 'YScale', 'log', 'XScale', 'log');
    else
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center'); axis off;
    end
    xticks([0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000]);
    xticklabels({'0.1','0.2','0.5','1','2','5','10','20','50','100','200','500','1k','2k','5k','10k','20k','50k','100k'});
    xlabel('ESD (µm)', 'FontSize', 8);
    ylabel('NBSS (mm^3/mm^3/m^3)', 'FontSize', 8);
    title('NBSS', 'FontSize', 10);

    % Panel 3 — relative abundance pie
    subplot(3,4,3);
    relAb_pie = relAb_group;
    relAb_pie(~isfinite(relAb_pie)) = 0;
    if sum(relAb_pie) > 0
        colormap(plankton_group_colors);
        pie(relAb_pie, repmat({''}, size(relAb_pie)));
    else
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center'); axis off;
    end
    title('Relative abundance', 'FontSize', 10);
    l = legend(plankton_groups_data_group, 'FontSize', 8, 'Location', 'west', 'NumColumns', 2);
    set(l, 'Position', [0.72 0.45 0.2 0.41]);
    title(l, 'Plankton groups', 'FontSize', 10);

    % Panel 5 — trophic pyramid
    subplot(3,4,5);
    color_troph = trophic_colors(plankton_groups_data_troph);
    height_t = str2double(plankton_groups_data_troph);
    width_t  = log(data_bv_troph(i,:) + 1);
    h_patch  = gobjects(numel(height_t), 1);
    for j = 1:numel(height_t)
        xp = [-width_t(j), width_t(j), width_t(j), -width_t(j)];
        yp = [height_t(j)-0.25, height_t(j)-0.25, height_t(j)+0.25, height_t(j)+0.25];
        h_patch(j) = patch(xp, yp, color_troph(j,:));
        hold on
    end
    xlabel('log biovolume +1 (mm^3/m^3)', 'FontSize', 8);
    title('Trophic pyramid', 'FontSize', 10);
    l = legend(h_patch, plankton_groups_data_troph_name, 'FontSize', 8, 'Location', 'west', 'NumColumns', 2);
    set(l, 'Position', [0.130 0.150 0.15 0.15]);
    title(l, 'Trophic groups', 'FontSize', 10);

    % Panel 6 — BSS
    ax6 = subplot(3,4,6);
    pos6 = get(ax6, 'Position');
    set(ax6, 'Position', [pos6(1)+0.03, pos6(2), pos6(3), pos6(4)]);
    bss_data = base(i).regroupped.Ybv_Ellipsoid_BV_spectra(:, base(i).regroupped.livingplace(1)) .* base(i).tot.X1;
    valid_bss = bss_data > 0;
    if any(valid_bss)
        plot(ESD(valid_bss), bss_data(valid_bss), ...
             'o', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
        set(gca, 'YScale', 'log', 'XScale', 'log');
    else
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center'); axis off;
    end
    xticks([0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000]);
    xticklabels({'0.1','0.2','0.5','1','2','5','10','20','50','100','200','500','1k','2k','5k','10k','20k','50k','100k'});
    xlim_bss = xlim(gca);
    xlabel('ESD (µm)', 'FontSize', 8);
    ylabel('BSS (mm^3/m^3)', 'FontSize', 8);
    title('BSS', 'FontSize', 10);

    % Panel 7 — relative biovolume pie
    subplot(3,4,7);
    bv_group_row = data_bv_group(i,:);
    bv_group_tot = sum(bv_group_row, 'omitnan');
    if bv_group_tot > 0
        relBv_pie = 100 * bv_group_row ./ bv_group_tot;
    else
        relBv_pie = zeros(size(bv_group_row));
    end
    relBv_pie(~isfinite(relBv_pie)) = 0;
    if sum(relBv_pie) > 0
        colormap(plankton_group_colors);
        pie(relBv_pie, repmat({''}, size(relBv_pie)));
    else
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center'); axis off;
    end
    title('Relative biovolume', 'FontSize', 10);

    % Panel 10 — biovolume composition per size class
    ax10 = subplot(3,4,10);
    pos10 = get(ax10, 'Position');
    set(ax10, 'Position', [pos10(1)+0.03, pos10(2), pos10(3), pos10(4)]);

    pp_i = base(i).regroupped.planktongroupplace;
    Ybv_denorm_group = base(i).regroupped.Ybv_Ellipsoid_BV_spectra(:, pp_i(1):pp_i(2)) .* base(i).tot.X1;
    Ybv_denorm_group = Ybv_denorm_group(:, col_order);
    % Apply same group exclusion as data_bv_group (excl_idx_group from filter_groups)
    if ~isempty(excl_idx_group)
        Ybv_denorm_group(:, excl_idx_group) = [];
    end
    Ybv_tot_persize = sum(Ybv_denorm_group, 2, 'omitnan');
    nonzero_rows    = Ybv_tot_persize > 0;
    if any(nonzero_rows)
        relative_Ybv = 100 * (Ybv_denorm_group(nonzero_rows,:) ./ Ybv_tot_persize(nonzero_rows));
        hbar_sz = bar(log10(ESD(nonzero_rows)), relative_Ybv, 'stacked', 'BarWidth', 1);
        for j = 1:size(plankton_group_colors, 1)
            if j <= numel(hbar_sz)
                hbar_sz(j).FaceColor = plankton_group_colors(j,:);
            end
        end
        set(hbar_sz, 'EdgeColor', 'none');
        ylim([0 100]);
        xlim(log10(xlim_bss));
        x_ticks = [0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000];
        set(gca, 'XTick', log10(x_ticks), 'XTickLabel', arrayfun(@num2str, x_ticks, 'UniformOutput', false));
    else
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center'); axis off;
    end
    xlabel('ESD (µm)', 'FontSize', 8);
    ylabel('Biovolume (mm^3/m^3)', 'FontSize', 8);
    title('Biovolume composition per size class', 'FontSize', 10);

    % Panel 11-12 — station map
    subplot(3,4,[11,12]);
    m_proj('mercator', 'lon', [lonmin lonmax], 'lat', [latmin latmax]);
    m_coast('patch', landcolor, 'edgecolor', [0.6 0.6 0.6]);
    m_grid('fontsize', 8, 'linestyle', 'none', 'meridianlabel', 'off', ...
           'parallellabel', 'off', 'xtick', [], 'ytick', []);
    hold on
    m_scatter(lon(i), lat(i), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.5 0.5 0.5]);

    if doSave
        outputDir = 'raw analysis per station';
        mkdir(outputDir);
        exportgraphics(gcf, fullfile(outputDir, [sid '.pdf']), ...
            'ContentType', 'vector', 'BackgroundColor', 'none');
    end
    pause(1); close gcf
end

%% ── 9. All-station summaries ─────────────────────────────────────────────

list3 = {'Abundance', 'Biovolume'};
[indx3, tf] = listdlg('PromptString', 'Data for bar charts and map per station', ...
    'ListString', list3, 'SelectionMode', 'single', 'ListSize', [400, 200]);
if ~tf, indx3 = 1; end

if indx3 == 1; data = data_ab_taxo; else; data = data_bv_taxo; end
n_samp = size(data, 1);

% ── H-bar: taxonomic level
p_group_colors;
num_groups    = numel(plankton_groups_data_taxo);
num_base_cols = size(plankton_group_colors_palet, 1);
cmap = zeros(num_groups, 3);
for i = 1:num_groups
    bi    = floor((i-1) * (num_base_cols-1) / max(num_groups-1, 1)) + 1;
    bi_nx = min(bi+1, num_base_cols);
    alpha = (i-1) * (num_base_cols-1) / max(num_groups-1, 1) - (bi-1);
    cmap(i,:) = (1-alpha) * plankton_group_colors_palet(bi,:) + alpha * plankton_group_colors_palet(bi_nx,:);
end

figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
if any(data(:))
    hbar = bar(data, 'stacked', 'FaceColor', 'flat');
    for k2 = 1:numel(hbar); hbar(k2).CData = k2; end
    colormap(cmap);
    set(gca, 'TickLength', [0.002 0.002], 'xtick', 1:n_samp, ...
        'xticklabel', sampleid, 'XTickLabelRotation', 75, ...
        'xlim', [0 n_samp+1], 'ylim', [0 max(sum(data,2))]);
    set(hbar, 'Edgecolor', [0 0 0]);
    l = legend(hbar, plankton_groups_data_taxo, 'FontSize', 8, 'Location', 'eastoutside', 'NumColumns', 2);
    title(l, 'Plankton groups', 'FontSize', 10);
end
if indx3 == 1; ylabel('Abundance (ind m^{-3})', 'FontSize', 10);
else;          ylabel('Biovolume (mm m^{-3})',   'FontSize', 10); end
if doSave
    mkdir('global raw analysis');
    if indx3==1; fname='Hbar taxo of abundance.pdf'; else; fname='Hbar taxo of biovolume.pdf'; end
    exportgraphics(gcf, fullfile('global raw analysis', fname), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(2); close gcf

% ── H-bar: PFT group level
if indx3 == 1; data = data_ab_group; else; data = data_bv_group; end
figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
if any(data(:))
    hbar = bar(data, 'stacked', 'FaceColor', 'flat');
    for k2 = 1:numel(hbar); hbar(k2).CData = k2; end
    colormap(plankton_group_colors);
    set(gca, 'TickLength', [0.002 0.002], 'xtick', 1:n_samp, ...
        'xticklabel', sampleid, 'XTickLabelRotation', 75, ...
        'xlim', [0 n_samp+1], 'ylim', [0 max(sum(data,2))]);
    set(hbar, 'Edgecolor', [0 0 0]);
    l = legend(hbar, plankton_groups_data_group, 'FontSize', 8, 'Location', 'eastoutside', 'NumColumns', 1);
    title(l, 'Plankton groups', 'FontSize', 10);
end
if indx3 == 1; ylabel('Abundance (ind m^{-3})', 'FontSize', 10);
else;          ylabel('Biovolume (mm m^{-3})',   'FontSize', 10); end
if doSave
    if indx3==1; fname='Hbar of abundance.pdf'; else; fname='Hbar of biovolume.pdf'; end
    exportgraphics(gcf, fullfile('global raw analysis', fname), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(2); close gcf

% ── H-bar: relative abundance/biovolume per PFT group
row_totals = sum(data, 2, 'omitnan');
row_totals(row_totals == 0) = NaN;
relative = 100 * bsxfun(@rdivide, data, row_totals);
relative(~isfinite(relative)) = 0;
figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
if any(relative(:))
    hbar = bar(relative, 'stacked', 'FaceColor', 'flat', 'BarWidth', 1);
    for k2 = 1:numel(hbar); hbar(k2).FaceColor = plankton_group_colors(k2,:); end
    set(gca, 'TickLength', [0.002 0.002], 'xtick', 1:n_samp, ...
        'xticklabel', strrep(sampleid, '_', '\_'), 'XTickLabelRotation', 75, 'ylim', [0 100]);
    set(hbar, 'Edgecolor', [0 0 0]);
    l = legend(hbar, plankton_groups_data_group, 'FontSize', 8, 'Location', 'eastoutside', 'NumColumns', 1);
    title(l, 'Plankton groups', 'FontSize', 10);
end
if doSave
    if indx3==1; fname='Relative Hbar of abundance.pdf'; else; fname='Relative Hbar of biovolume.pdf'; end
    exportgraphics(gcf, fullfile('global raw analysis', fname), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(2); close gcf

% ── Map: total abundance or biovolume
data_sum = sum(data, 2, 'omitnan');
figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
m_proj('mercator', 'lon', [lonmin lonmax], 'lat', [latmin latmax]);
m_coast('patch', landcolor, 'edgecolor', [0.6 0.6 0.6]); hold on
m_grid('fontsize', 8, 'linestyle', 'none', 'meridianlabel', 'off', 'parallellabel', 'off', 'xtick', [], 'ytick', []);
valid_map = data_sum > 0;
if any(valid_map)
    m_scatter(lon(valid_map), lat(valid_map), 50, log10(data_sum(valid_map)), ...
              'filled', 'markeredgecolor', [0.5 0.5 0.5]);
    load('slanCM_Data.mat');
    color = slandarerCM(6).Colors{1,9};
    c2 = colorbar; colormap(flipud(color));
    cur_ticks = c2.Ticks;
    new_labels = 10.^cur_ticks;
    is_pow10 = abs(log10(new_labels) - round(log10(new_labels))) < 1e-10;
    if any(is_pow10)
        set(c2, 'Ticks', cur_ticks(is_pow10), ...
            'TickLabels', arrayfun(@num2str, new_labels(is_pow10), 'UniformOutput', false));
    end
    if indx3==1; title(c2, 'Abundance (ind m^{-3})', 'FontSize', 11, 'FontWeight', 'bold');
    else;        title(c2, 'Biovolume (mm m^{-3})',   'FontSize', 11, 'FontWeight', 'bold'); end
end
if doSave
    mkdir('global raw analysis');
    if indx3==1; fname='Map of abundance.pdf'; else; fname='Map of biovolume.pdf'; end
    exportgraphics(gcf, fullfile('global raw analysis', fname), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(2); close gcf

% ── Map: Shannon diversity
load('slanCM_Data.mat');
diversityH = zeros(m, 1);
for i = 1:m
    ab_row = data_ab_taxo(i,:);
    tot_i  = sum(ab_row, 'omitnan');
    if tot_i > 0
        rel = ab_row ./ tot_i;
        diversityH(i) = -sum(rel .* log(rel), 'omitnan');
    end
end
figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
m_proj('mercator', 'lon', [lonmin lonmax], 'lat', [latmin latmax]);
m_coast('patch', landcolor, 'edgecolor', [0.6 0.6 0.6]); hold on
m_grid('fontsize', 8, 'linestyle', 'none', 'meridianlabel', 'off', 'parallellabel', 'off', 'xtick', [], 'ytick', []);
color2 = slandarerCM(2).Colors{1,15};
colormap(flipud(color2));
m_scatter(lon, lat, 50, diversityH, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
c = colorbar; title(c, 'Shannon Index', 'FontSize', 11, 'FontWeight', 'bold');
if doSave
    exportgraphics(gcf, fullfile('global raw analysis', 'Map of diversity.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(2); close gcf

%% ── 10. Optional per-group analysis ──────────────────────────────────────

list5 = {'Yes', 'No'};
[indx5, tf] = listdlg('PromptString', 'Select one taxon/group for dedicated analysis?', ...
    'ListString', list5, 'SelectionMode', 'single', 'ListSize', [400, 200]);

if tf && indx5 == 1

    list3b = {'Abundance', 'Biovolume'};
    [indx3b, tf] = listdlg('PromptString', 'Data for H-Bar and Map for one group', ...
        'ListString', list3b, 'SelectionMode', 'single', 'ListSize', [400, 200]);
    if ~tf, indx3b = 1; end

    list1 = {'one taxon'; 'one plankton group'; 'one trophic group'; ...
             'all taxon individually'; 'all plankton group'; 'all trophic group'};
    [indx1, tf] = listdlg('PromptString', 'Choose grouping level', ...
        'ListString', list1, 'SelectionMode', 'single', 'ListSize', [400, 200]);
    if ~tf, indx1 = 1; end

    switch indx1
        case 1
            list2b = plankton_groups_data_taxo;
            [indx2, tf] = listdlg('PromptString', 'Choose taxon', 'ListString', list2b, 'ListSize', [400, 200]);
            if ~tf, indx2 = 1; end
            plankton_group = list2b(indx2);
            if indx3b==1; data_one = data_ab_taxo(:,indx2); else; data_one = data_bv_taxo(:,indx2); end
        case 2
            list2b = plankton_groups_data_group;
            [indx2, tf] = listdlg('PromptString', 'Choose group', 'ListString', list2b, 'ListSize', [400, 200]);
            if ~tf, indx2 = 1; end
            plankton_group = list2b(indx2);
            if indx3b==1; data_one = data_ab_group(:,indx2); else; data_one = data_bv_group(:,indx2); end
        case 3
            list2b = plankton_groups_data_troph_name;
            [indx2, tf] = listdlg('PromptString', 'Choose trophic group', 'ListString', list2b, 'ListSize', [400, 200]);
            if ~tf, indx2 = 1; end
            plankton_group = list2b(indx2);
            if indx3b==1; data_one = data_ab_troph(:,indx2); else; data_one = data_bv_troph(:,indx2); end
        case 4
            plankton_group = plankton_groups_data_taxo;
            if indx3b==1; data_one = data_ab_taxo; else; data_one = data_bv_taxo; end
        case 5
            plankton_group = plankton_groups_data_group;
            if indx3b==1; data_one = data_ab_group; else; data_one = data_bv_group; end
        case 6
            plankton_group = plankton_groups_data_troph_name;
            if indx3b==1; data_one = data_ab_troph; else; data_one = data_bv_troph; end
    end

    plankton_group = cellstr(plankton_group);

    for gi = 1:size(data_one, 2)

        figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
        col_data = data_one(:, gi);
        if any(col_data)
            hbar_one = bar(col_data, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0 0 0]);
        end
        set(gca, 'TickLength', [0.002 0.002], 'xtick', 1:n_samp, ...
            'xticklabel', sampleid, 'XTickLabelRotation', 75, ...
            'xlim', [0 n_samp+1], 'ylim', [0 max(col_data)]);
        if indx3b==1; ylabel('Abundance (ind m^{-3})', 'FontSize', 10);
        else;          ylabel('Biovolume (mm m^{-3})',  'FontSize', 10); end
        legend(plankton_group(gi), 'FontSize', 8, 'Location', 'eastoutside');

        if doSave
            if indx3b==1
                outputDir = fullfile('global raw analysis', ['Hbar abundance for ' list1{indx1}]);
            else
                outputDir = fullfile('global raw analysis', ['Hbar biovolume for ' list1{indx1}]);
            end
            mkdir(outputDir);
            exportgraphics(gcf, fullfile(outputDir, ['Hbar ' char(plankton_group{gi}) '.pdf']), ...
                'ContentType', 'vector', 'BackgroundColor', 'none');
        end
        pause(1); close gcf

        % Map for this group
        datamap = col_data;
        figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
        m_proj('mercator', 'lon', [lonmin lonmax], 'lat', [latmin latmax]);
        m_coast('patch', landcolor, 'edgecolor', [0.6 0.6 0.6]); hold on
        m_grid('fontsize', 8, 'linestyle', 'none', 'meridianlabel', 'off', 'parallellabel', 'off', 'xtick', [], 'ytick', []);
        if max(datamap) > min(datamap)
            scaled_size = 10 + (datamap - min(datamap)) / (max(datamap) - min(datamap)) * 90;
        else
            scaled_size = repmat(50, size(datamap));
        end
        m_scatter(lon, lat, scaled_size, datamap, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
        c2 = colorbar;
        if indx3b==1; title(c2, 'Abundance (ind m^{-3})', 'FontSize', 10, 'FontWeight', 'bold');
        else;          title(c2, 'Biovolume (mm m^{-3})',  'FontSize', 10, 'FontWeight', 'bold'); end

        if doSave
            if indx3b==1
                outputDir = fullfile('global raw analysis', ['Map abundance for ' list1{indx1}]);
            else
                outputDir = fullfile('global raw analysis', ['Map biovolume for ' list1{indx1}]);
            end
            mkdir(outputDir);
            exportgraphics(gcf, fullfile(outputDir, ['Map ' char(plankton_group{gi}) '.pdf']), ...
                'ContentType', 'vector', 'BackgroundColor', 'none');
        end
        pause(1); close gcf
    end
end

%% ── 11. NBSS analysis ────────────────────────────────────────────────────

mkdir(fullfile(pwd, 'NBSS analysis'));

% All-station NBSS
figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
for i = 1:numel(base)
    d = base(i).regroupped.NBSSliving;
    valid = d > 0;
    if any(valid)
        hold on
        plot(ESD(valid), d(valid), 'o', ...
            'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    end
end
set(gca, 'YScale', 'log', 'XScale', 'log');
xticks([0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'0.1','0.2','0.5','1','2','5','10','20','50','100','200','500','1k','2k','5k','10k','20k','50k','100k'});
xlabel('ESD (µm)', 'FontSize', 10); ylabel('NBSS (mm^3/mm^3/m^3)', 'FontSize', 10);
title('NBSS — all stations', 'FontSize', 10);
if doSave
    exportgraphics(gcf, fullfile('NBSS analysis', 'NBSS Spectra.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1); close gcf

% Truncated NBSS
figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
for i = 1:numel(base)
    d = base(i).regroupped.NBSSliving_trunked;
    valid = d > 0 & isfinite(d);
    if any(valid)
        hold on
        plot(ESD(valid), d(valid), 'o', ...
            'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    end
end
set(gca, 'YScale', 'log', 'XScale', 'log');
xticks([0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'0.1','0.2','0.5','1','2','5','10','20','50','100','200','500','1k','2k','5k','10k','20k','50k','100k'});
xlabel('ESD (µm)', 'FontSize', 10); ylabel('NBSS (mm^3/mm^3/m^3)', 'FontSize', 10);
title('NBSS truncated — all stations', 'FontSize', 10);
if doSave
    exportgraphics(gcf, fullfile('NBSS analysis', 'NBSS Spectra trunked.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1); close gcf

% Map NBSS slopes
slope_values = arrayfun(@(b) b.regroupped.NBSSliving_slope, base);
figure('NumberTitle', 'off', 'Position', [10 50 scrsz(3) scrsz(4)-150], 'Color', 'white');
m_proj('mercator', 'lon', [lonmin lonmax], 'lat', [latmin latmax]);
m_coast('patch', landcolor, 'edgecolor', [0.6 0.6 0.6]); hold on
m_grid('fontsize', 8, 'linestyle', 'none', 'meridianlabel', 'off', 'parallellabel', 'off', 'xtick', [], 'ytick', []);
m_scatter(lon, lat, 50, slope_values, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
load('slanCM_Data.mat');
color = slandarerCM(6).Colors{1,9};
colormap(flipud(color));
c = colorbar; title(c, 'NBSS Slope', 'FontWeight', 'bold');
if doSave
    exportgraphics(gcf, fullfile('NBSS analysis', 'Map of NBSS Slope.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1); close gcf

%% ── Local helper functions ───────────────────────────────────────────────

% filter_groups: interactive exclusion dialog
% Returns also excl_idx: indices of excluded columns (in original order before exclusion)
function [data_ab, data_bv, groups, excl_idx] = filter_groups(data_ab, data_bv, groups, label)
    excl_idx = [];
    list2 = {'Exclude manually', 'Load exclusion list', 'No exclusion'};
    [indx, tf2] = listdlg('PromptString', ['Exclude some ' label ' from analysis?'], ...
        'ListString', list2, 'SelectionMode', 'single', 'ListSize', [400, 200]);
    if ~tf2, return; end
    switch indx
        case 1
            [indx2, tf3] = listdlg('PromptString', 'Select groups to exclude (Ctrl/Cmd for multiple)', ...
                'ListString', groups, 'ListSize', [400, 200]);
            if ~tf3, return; end
            to_excl = num2cell(indx2(:));
            to_excl(:,2) = groups(indx2);
            writetable(cell2table(to_excl), [label '_excluded.xlsx'], 'WriteVariableNames', false);
            excl_idx = indx2(:);
            data_ab(:, indx2) = [];
            data_bv(:, indx2) = [];
            groups(indx2)     = [];
        case 2
            f = msgbox(['Select the exclusion list for ' label]);
            uiwait(f);
            [file, path] = uigetfile('*.xlsx');
            T = readtable(fullfile(path, file), 'ReadVariableNames', false);
            indx2 = table2array(T(:,1));
            excl_idx = indx2(:);
            data_ab(:, indx2) = [];
            data_bv(:, indx2) = [];
            groups(indx2)     = [];
    end
end

% trophic_colors: returns RGB matrix for a trophic label cell array
function cmat = trophic_colors(g)
    lut = {'-1',[1 1 1]; '1',[0.40 0.69 0.39]; '1.5',[0.63 0.77 0.53]; ...
           '2',[0.3 0.54 0.69]; '2.5',[1 0.93 0.66]; '3',[0.78 0.33 0.149]; '3.5',[0.9 0.9 0.9]};
    cmat = repmat([0.8 0.8 0.8], numel(g), 1);
    for gi = 1:numel(g)
        row = find(strcmp(strtrim(g{gi}), lut(:,1)), 1);
        if ~isempty(row), cmat(gi,:) = lut{row,2}; end
    end
end
