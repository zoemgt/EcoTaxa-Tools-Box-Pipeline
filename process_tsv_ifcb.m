% process_tsv_ifcb.m
% Assembles data from EcoTaxa TSV files (IFCB) into the structured "base" database.

function [base, m, zoo_groups, Idlistshort] = process_tsv_ifcb

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
    base(i).sampleid  = sample;
    base(i).date      = unique(S.object_date);
    base(i).time      = unique(S.object_time);
    base(i).latitude  = double(unique(S.object_lat));
    base(i).longitude = double(unique(S.object_lon));

    if ismember('object_lat_end', S.Properties.VariableNames)
        base(i).latitude_end = unique(S.object_lat_end);
    end
    if ismember('object_lon_end', S.Properties.VariableNames)
        base(i).longitude_end = unique(S.object_lon_end);
    end

    base(i).vol = toNum(unique(S.acq_volume_sampled)) / 1e6;  % mL -> m3

    %% Pixel size [mm/pixel]
    pixelsize = (1 / toNum(unique(S.acq_resolution_pixel_per_micron))) / 1000;
    base(i).pixelsize = pixelsize;

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
        base(i).(fn).major            = S.object_major_axis_length(I)  * pixelsize;
        base(i).(fn).minor            = S.object_minor_axis_length(I)  * pixelsize;
        base(i).(fn).area             = S.object_surface_area(I)       * pixelsize^2;
        base(i).(fn).ESD              = 2 * sqrt((S.object_surface_area(I) * pixelsize^2) / pi);
        base(i).(fn).summedbiovolume  = S.object_summed_biovolume(I)         * pixelsize^3;
        base(i).(fn).summedarea       = S.object_summed_surface_area(I)      * pixelsize^2;

        % Conversion factor [m-3]
        base(i).(fn).conver = 1 ./ base(i).vol;

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
