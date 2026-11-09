% process_tsv_planktoscope.m
% Assembles data from EcoTaxa TSV files (PlanktoScope) into the structured "base" database.

function [base, m, zoo_groups, Idlistshort] = process_tsv_planktoscope

%% Initialise

A         = dir('*.tsv');
filenames = {A.name};
m         = numel(filenames);

base        = [];
Idlist      = {};
Idlistshort = {};

h = waitbar(0, 'Loading TSV files and building base...');

for i = 1:m

    waitbar(i/m, h, sprintf('Loading TSV files: %d of %d', i, m));

    S      = readtable(filenames{i}, 'FileType', 'text', 'ReadVariableNames', true);
    sample = unique(S.sample_id);

    %% Homogenise taxonomic annotation hierarchy (replace '-' and '>' with '_')
    S.object_annotation_hierarchy = strrep(S.object_annotation_hierarchy, '-', '_');
    S.object_annotation_hierarchy = strrep(S.object_annotation_hierarchy, '>', '_');

    %% Build base — metadata
    base(i).sampleid       = sample;
    base(i).ship           = unique(S.sample_ship);
    base(i).scientificprog = unique(S.sample_project);
    base(i).date           = unique(S.object_date);
    base(i).time           = unique(S.object_time);
    base(i).latitude       = double(unique(S.object_lat));
    base(i).longitude      = double(unique(S.object_lon));

    if ismember('object_lat_end', S.Properties.VariableNames)
        base(i).latitude_end = unique(S.object_lat_end);
    end
    if ismember('object_lon_end', S.Properties.VariableNames)
        base(i).longitude_end = unique(S.object_lon_end);
    end

    base(i).nettype = unique(S.sample_sampling_gear);
    base(i).zmax    = unique(S.object_depth_max);
    base(i).zmin    = unique(S.object_depth_min);

    %% Pixel size [mm/pixel]
    rawpix = unique(S.process_pixel);
    if isnumeric(rawpix)
        pixelsize = rawpix / 1000;
    else
        pixelsize = str2double(cell2mat(rawpix)) / 1000;
    end
    base(i).pixelsize = pixelsize;

    %% Volume-related fields
    % readtable in 2025b may import numeric columns as cell arrays of strings
    % when the TSV contains mixed or missing values — toNum() handles both cases.
    base(i).volconc   = toNum(unique(S.sample_concentrated_sample_volume)) / 1e6; % mL -> m3
    base(i).nb_frame  = toNum(unique(S.acq_nb_frame));
    base(i).celltype  = toNum(unique(S.acq_celltype)) / 1e6;                       % µm -> m
    base(i).resolution = unique(S.acq_camera_resolution);

    base(i).sample_initial_col_vol_m3 = toNum(unique(S.sample_total_volume)) / 1000; % L -> m3
    base(i).volume_imaged_processed   = toNum(unique(S.acq_imaged_volume))   / 1e6;  % mL -> m3

    % Dilution / concentration correction (multiplicative factor applied to conver)
    % Priority 1: sample_dilution_factor (multiplied)
    % Priority 2: sample_concentration_factor (divided)
    % Default: 1 (no correction)
    base(i).dilution_factor      = NaN;
    base(i).concentration_factor = NaN;
    if ismember('sample_dilution_factor', S.Properties.VariableNames)
        df = toNum(unique(S.sample_dilution_factor));
        if ~isempty(df) && ~all(isnan(df)), base(i).dilution_factor = df(1); end
    end
    if ismember('sample_concentration_factor', S.Properties.VariableNames)
        cf = toNum(unique(S.sample_concentration_factor));
        if ~isempty(cf) && ~all(isnan(cf)), base(i).concentration_factor = cf(1); end
    end

    if ~isnan(base(i).dilution_factor)
        base(i).volume_factor = base(i).dilution_factor;           % multiply by dilution factor
    elseif ~isnan(base(i).concentration_factor)
        base(i).volume_factor = 1 / base(i).concentration_factor;  % divide by concentration factor
    else
        base(i).volume_factor = 1;                                 % no correction
    end

    %% Size fractions
    base(i).fracids = unique(S.acq_id);
    nfrac           = numel(base(i).fracids);

    for fracnb = 1:nfrac

        fn = ['d' num2str(fracnb)];
        I  = strcmp(S.acq_id, base(i).fracids{fracnb});

        % Taxonomic annotation
        base(i).(fn).object_annotation_hierarchy = S.object_annotation_hierarchy(I);

        % Update global taxonomic lists
        [Idlist, ia, ~] = unique([Idlist; S.object_annotation_hierarchy(I)]);
        Idlistshort     = [Idlistshort; S.object_annotation_category(I)];
        Idlistshort     = Idlistshort(ia);

        % Percent validated
        isVal = ~cellfun(@isempty, strfind(S.object_annotation_status(I), 'validated'));
        base(i).(fn).percentValidated = 100 * sum(isVal) / numel(isVal);

        % Morphological measurements [mm or mm2]
        if isnumeric(S.object_major)
            base(i).(fn).major    = S.object_major(I)    * pixelsize;
            base(i).(fn).minor    = S.object_minor(I)    * pixelsize;
            base(i).(fn).area_exc = S.object_area_exc(I) * pixelsize^2;
            base(i).(fn).area     = S.object_area(I)     * pixelsize^2;
            base(i).(fn).ESD      = 2 * sqrt((S.object_area(I) * pixelsize^2) / pi);
        else
            base(i).(fn).major    = cellfun(@str2double, S.object_major(I))    * pixelsize;
            base(i).(fn).minor    = cellfun(@str2double, S.object_minor(I))    * pixelsize;
            base(i).(fn).area_exc = cellfun(@str2double, S.object_area_exc(I)) * pixelsize^2;
            base(i).(fn).area     = cellfun(@str2double, S.object_area(I))     * pixelsize^2;
            base(i).(fn).ESD      = 2 * sqrt((cellfun(@str2double, S.object_area(I)) * pixelsize^2) / pi);
        end

        % Conversion factor [m-3]
        base(i).(fn).conver = base(i).volconc .* base(i).volume_factor ./ ...
            (base(i).volume_imaged_processed .* base(i).sample_initial_col_vol_m3);

    end
end

close(h);

m          = numel(base);
zoo_groups = Idlist;

end

%% Helper — convert numeric or cell-of-strings scalar to double
function v = toNum(x)
    if isnumeric(x)
        v = x;
    elseif iscell(x)
        v = str2double(cell2mat(x));
    else
        v = str2double(x);
    end
end
