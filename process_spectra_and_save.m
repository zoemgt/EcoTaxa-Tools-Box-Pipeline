% process_spectra_and_save.m
% Computes abundance, biovolume, and NBSS spectra from the "base" struct,
% assembles multi-fraction samples, applies taxonomic grouping,
% computes NBSS slopes, and saves the final .mat base.

%% Spectral computation settings
smin = 1e-12;   % lower biovolume limit
smax = 10000;   % upper biovolume limit (ESD ~26 mm)
k    = 2^(1/4); % logarithmic bin base (Platt & Denman, 1978)
uu   = 1;       % first size class for regression

%% Compute abundance, biovolume and NBSS for each fraction
h = waitbar(0, 'Computing abundance, biovolume and NBSS spectra...');

for i = 1:m
    nfrac = numel(base(i).fracids);
    waitbar(i/m, h, sprintf('Computing spectra: %d of %d', i, m));

    for fracnb = 1:nfrac
        subset = base(i).(['d' num2str(fracnb)]);
        [base, ~] = process_abundances_biovolumes_NBSSspectra( ...
            subset, base, smin, smax, k, uu, zoo_groups, i, fracnb);
    end
end
close(h);

%% Assemble fractions into the "tot" level
for i = 1:m
    nfrac = numel(base(i).fracids);

    Ab                          = [];
    Bv                          = [];
    Ybv_Plain_Area_BV_spectra   = [];
    Ybv_Riddled_Area_BV_spectra = [];
    Ybv_Ellipsoid_BV_spectra    = [];
    Yab                         = [];

    for fracnb = 1:nfrac
        fn = ['d' num2str(fracnb)];
        d  = base(i).(fn);

        if fracnb == 1
            Ab                          = d.Ab;
            Bv                          = d.Bv;
            Ybv_Plain_Area_BV_spectra   = d.Ybv_Plain_Area_BV_spectra;
            Ybv_Riddled_Area_BV_spectra = d.Ybv_Riddled_Area_BV_spectra;
            Ybv_Ellipsoid_BV_spectra    = d.Ybv_Ellipsoid_BV_spectra;
            Yab                         = d.Yab;
        else
            Ab  = sum([Ab,  d.Ab],  2, 'omitnan');
            Bv  = sum([Bv,  d.Bv],  2, 'omitnan');
            Ybv_Plain_Area_BV_spectra   = sum(cat(3, Ybv_Plain_Area_BV_spectra,   d.Ybv_Plain_Area_BV_spectra),   3, 'omitnan');
            Ybv_Riddled_Area_BV_spectra = sum(cat(3, Ybv_Riddled_Area_BV_spectra, d.Ybv_Riddled_Area_BV_spectra), 3, 'omitnan');
            Ybv_Ellipsoid_BV_spectra    = sum(cat(3, Ybv_Ellipsoid_BV_spectra,    d.Ybv_Ellipsoid_BV_spectra),    3, 'omitnan');
            Yab = sum(cat(3, Yab, d.Yab), 3, 'omitnan');
        end
    end

    base(i).tot.Ab                          = Ab;
    base(i).tot.Bv                          = Bv;
    base(i).tot.Ybv_Plain_Area_BV_spectra   = Ybv_Plain_Area_BV_spectra;
    base(i).tot.Ybv_Riddled_Area_BV_spectra = Ybv_Riddled_Area_BV_spectra;
    base(i).tot.Ybv_Ellipsoid_BV_spectra    = Ybv_Ellipsoid_BV_spectra;
    base(i).tot.Yab                         = Yab;

    % Copy size-class vectors from last fraction
    lastfn = ['d' num2str(fracnb)];
    base(i).tot.X          = base(i).(lastfn).X;
    base(i).tot.X1         = base(i).(lastfn).X1;
    base(i).tot.Zoo_groups = base(i).(lastfn).Zoo_groups;

    % ESD vector [µm]
    sizelist = exp(base(i).tot.X);          % mm3
    ESD      = 2 * (sizelist * 3 / (4*pi)).^(1/3) * 1000;  % mm -> µm
    base(i).tot.ESDvector = ESD;
end

%% Load trophic affiliation table
table_groupage = readtable('trophic_affiliation_of_organisms.xlsx', 'ReadVariableNames', false);
table_groupage = table2cell(table_groupage);

msgbox({'-1  do not feed'; '1   phototrophs'; '1.5 mixotrophs'; ...
        '2   grazers'; '2.5 omnivorous'; '3   predators'; '3.5 unknown trophic group'}, ...
        'Trophic group legend');

Zoo_groups_ref = table_groupage(:,1);

%% Handle temporary groups (t00X)
istemporary = false;
k_temp = cellfun(@(x) ~isempty(strfind(x, 'temporary_')), zoo_groups);

if any(k_temp)

    answer = questdlg( ...
        'Your files include one or several temporary "t00X" categories. Do you have a functional/trophic mapping for them?', ...
        'Temporary categories mapping', ...
        'Yes — load them', 'No — create them', 'No — ignore them (not recommended)', ...
        'Yes — load them');

    switch answer
        case 'No — ignore them (not recommended)'
            zoo_groups(k_temp) = [];
            for i = 1:m
                base(i).tot.Zoo_groups(k_temp)                      = [];
                base(i).tot.Ab(k_temp)                              = [];
                base(i).tot.Bv(k_temp)                              = [];
                base(i).tot.Ybv_Plain_Area_BV_spectra(:, k_temp)   = [];
                base(i).tot.Ybv_Riddled_Area_BV_spectra(:, k_temp) = [];
                base(i).tot.Ybv_Ellipsoid_BV_spectra(:, k_temp)    = [];
                base(i).tot.Yab(:, k_temp)                         = [];
            end

        case 'Yes — load them'
            f = msgbox('Select the temporary mapping *.xlsx');
            uiwait(f);
            [file, path] = uigetfile('*.xlsx');
            addontemp = table2cell(readtable(fullfile(path, file), 'ReadVariableNames', false));
            Zoo_groups_ref = [Zoo_groups_ref; addontemp(:,1)];
            istemporary = true;

        case 'No — create them'
            groups    = table_groupage(:, 2:end);
            group1    = unique(groups(:,1));
            group2    = unique(groups(:,2));
            group3    = unique(cellstr(num2str(cell2mat(groups(:,3)))));
            group1    = [{'living'};    setdiff(group1, 'living',       'stable')];
            group2    = [{'unidentified'}; setdiff(group2, 'unidentified', 'stable')];
            group3    = [{'3.5'};       setdiff(group3, '3.5',          'stable')];
            new_taxa  = zoo_groups(k_temp);
            newfunctional = [];
            for i = 1:numel(new_taxa)
                settings = settingsdlg( ...
                    'Description', ['New temporary category found: ' char(new_taxa(i))], ...
                    'title', 'New taxa functional mapping', ...
                    'Alive', group1, 'plankton group', group2, 'trophic group', group3, ...
                    'WindowWidth', 800);
                newfunctional = [newfunctional, settings]; %#ok<AGROW>
            end
            newfunctional = struct2table(newfunctional);
            addontemp = [new_taxa, table2cell(newfunctional)];
            tosave = array2table(addontemp);
            [file, path] = uiputfile('temporarymapping_instrument_net_location.xlsx');
            writetable(tosave, fullfile(path, file), 'WriteVariableNames', false);
            Zoo_groups_ref = [Zoo_groups_ref; addontemp(:,1)];
            istemporary = true;
    end
end

%% Handle new taxa not yet in the affiliation table
new_taxa = {};
for i = 1:numel(zoo_groups)
    if ~any(strcmp(zoo_groups{i}, Zoo_groups_ref))
        new_taxa{end+1,1} = zoo_groups{i}; %#ok<AGROW>
    end
end

if ~isempty(new_taxa)

    f = warndlg({['New taxa were found (' num2str(numel(new_taxa)) '). A window will open for each of them.'], ...
        'Your choices will be saved permanently in trophic_affiliation_of_organisms.xlsx and applied to all future runs.', ...
        'A backup copy of the current file will be created.'}, 'New taxa');
    uiwait(f);

    groups = table_groupage(:, 2:end);
    group1 = unique(groups(:,1));
    group2 = unique(groups(:,2));
    group3 = unique(cellstr(num2str(cell2mat(groups(:,3)))));
    newfunctional = [];
    for i = 1:numel(new_taxa)
        settings = settingsdlg( ...
            'Description', ['New taxonomic group found: ' char(new_taxa(i))], ...
            'title', 'New taxa functional mapping', ...
            'Alive', group1, 'plankton group', group2, 'trophic group', group3, ...
            'WindowWidth', 800);
        newfunctional = [newfunctional, settings]; %#ok<AGROW>
    end
    newfunctional = struct2table(newfunctional);
    addon = [new_taxa, table2cell(newfunctional)];
    if size(addon,2) >= 5, addon(:,5) = []; end
    table_groupage = [table_groupage; addon];
    tosave = array2table(table_groupage);

    copyfile('trophic_affiliation_of_organisms.xlsx', ...
        ['trophic_affiliation_of_organisms_backup_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.xlsx']);

    cd(directoryoftoolbox);
    writetable(tosave, 'trophic_affiliation_of_organisms.xlsx', 'WriteVariableNames', false);
    cd(folder);
end

if istemporary && size(addontemp,2) >= 5
    addontemp(:,5) = [];
end
if istemporary
    table_groupage = [table_groupage; addontemp];
end

%% Regroup taxa by living, PFT, and trophic levels
for i = 1:m
    base(i).regroupped = f_regroup_all(table_groupage, base(i).tot);
end

%% Summary resume tables
Ab_resume  = [];
Bv_resume  = [];
samplelist = {};

for i = 1:m
    Ab_resume  = [Ab_resume;  base(i).regroupped.Ab];
    Bv_resume  = [Bv_resume;  base(i).regroupped.Bv];
    samplelist = [samplelist; base(i).sampleid];
    base(i).regroupped.Idlistshort = Idlistshort;
end

%% NBSS extraction and truncation
% Active method: Ellipsoid (recommended, Vandromme et al. 2012; change field name below to switch method)
for i = 1:numel(base)
    r = base(i).regroupped;
    base(i).regroupped.NBSSliving  = r.Ybv_Ellipsoid_BV_spectra(:, r.livingplace);
    base(i).regroupped.NBSSpft     = r.Ybv_Ellipsoid_BV_spectra(:, r.planktongroupplace(1):r.planktongroupplace(2));
    base(i).regroupped.NBSStrophic = r.Ybv_Ellipsoid_BV_spectra(:, r.trophicplace(1):r.trophicplace(2));
end

% Helper: truncate NBSS (remove rising limb, zeros, isolated tails)
    function data_out = truncate_nbss(data_in)
        data_out = data_in;
        for col = 1:size(data_in, 2)
            d = data_in(:, col);
            [~, I] = max(d, [], 'omitnan');
            d(1:I)       = NaN;
            d(d == 0)    = NaN;
            I2 = find(diff([0; isnan(d)]) == 1);
            nanlength = zeros(numel(I2), 1);
            for j = 2:(numel(I2)-1)
                nanlength(j) = sum(isnan(d(I2(j):I2(j+1)-1)));
            end
            I3 = nanlength > 5;
            if any(I3)
                d(I2(find(I3,1,'first')):end) = NaN;
            end
            data_out(:, col) = d;
        end
    end

for i = 1:numel(base)
    base(i).regroupped.NBSSliving_trunked  = truncate_nbss(base(i).regroupped.NBSSliving);
    base(i).regroupped.NBSSpft_trunked     = truncate_nbss(base(i).regroupped.NBSSpft);
    base(i).regroupped.NBSStrophic_trunked = truncate_nbss(base(i).regroupped.NBSStrophic);
end

%% NBSS slope and intercept (linear regression of log10(NBSS) on log10(biovolume))
for i = 1:numel(base)
    X    = base(i).tot.X / log(10);   % class midpoints: ln(mm3) -> log10(mm3), same base as log10(NBSS)
    NBSS = log10(base(i).regroupped.NBSSliving_trunked);
    valid = ~isnan(X) & ~isnan(NBSS);

    if sum(valid) >= 5
        coef = polyfit(X(valid), NBSS(valid), 1);
        base(i).regroupped.NBSSliving_slope     = coef(1);   % slope b
        base(i).regroupped.NBSSliving_intercept = coef(2);   % intercept a = log10(NBSS) at 1 mm3
    else
        base(i).regroupped.NBSSliving_slope     = NaN;
        base(i).regroupped.NBSSliving_intercept = NaN;
    end
end

%% Save
prompt   = {'Instrument:', 'Project:'};
answer   = inputdlg(prompt, 'Save base as:', 1, {'instrument', 'project'});
savename = ['base_' answer{1} '_' answer{2}];
save(savename, 'base');
