% process_tsv_zooscan.m
% Assembles data from EcoTaxa TSV files (ZooScan) into the structured "base" database.

function [base, m, zoo_groups, Idlistshort] = process_tsv_zooscan

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
    base(i).scanop         = unique(S.sample_scan_operator);
    base(i).ship           = unique(S.sample_ship);
    base(i).scientificprog = unique(S.sample_program);
    base(i).stationid      = unique(S.sample_stationid);
    base(i).date           = unique(S.object_date);
    base(i).time           = unique(S.object_time);
    base(i).latitude       = double(unique(S.object_lat));
    base(i).longitude      = double(unique(S.object_lon));
    base(i).depth          = unique(S.sample_bottomdepth);
    base(i).CTDref         = unique(S.sample_ctdrosettefilename);
    base(i).otherref       = unique(S.sample_other_ref);
    base(i).barcode        = unique(S.sample_barcode);
    base(i).nettype        = unique(S.sample_net_type);
    base(i).townb          = unique(S.sample_tow_nb);
    base(i).towtype        = unique(S.sample_tow_type);
    base(i).netmesh        = unique(S.sample_net_mesh);
    base(i).netsurf        = unique(S.sample_net_surf);
    base(i).zmax           = unique(S.sample_zmax);
    base(i).zmin           = unique(S.sample_zmin);
    base(i).vol            = toNum(unique(S.sample_tot_vol));  % volume filtered [m3]
    base(i).samplecomments = unique(S.sample_comment);

    %% Pixel size [mm/pixel]
    pixelsize      = toNum(unique(S.process_particle_pixel_size_mm));
    base(i).pixelsize = pixelsize;

    %% Size fractions (d1 = >1000 µm, d2 = 100-1000 µm, or tot for no fractionation)
    base(i).fracids = unique(S.acq_id);
    nfrac           = numel(base(i).fracids);

    for fracnb = 1:nfrac

        fn = ['d' num2str(fracnb)];  % dynamic field name — replaces eval
        I  = strcmp(S.acq_id, base(i).fracids{fracnb});

        base(i).(fn).Fracmin    = unique(S.acq_min_mesh(I));
        base(i).(fn).Fracsup    = unique(S.acq_max_mesh(I));
        base(i).(fn).Fracnb     = unique(S.acq_sub_part(I));
        base(i).(fn).Resolution = unique(S.process_img_resolution(I));

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
            base(i).(fn).major     = S.object_major(I)     * pixelsize;
            base(i).(fn).minor     = S.object_minor(I)     * pixelsize;
            base(i).(fn).area_exc  = S.object_area_exc(I)  * pixelsize^2;
            base(i).(fn).area      = S.object_area(I)      * pixelsize^2;
            base(i).(fn).perimferet = S.object_feret(I)    * pixelsize;
            base(i).(fn).ESD       = 2 * sqrt((S.object_area(I) * pixelsize^2) / pi);
        else
            base(i).(fn).major     = cellfun(@str2double, S.object_major(I))     * pixelsize;
            base(i).(fn).minor     = cellfun(@str2double, S.object_minor(I))     * pixelsize;
            base(i).(fn).area_exc  = cellfun(@str2double, S.object_area_exc(I))  * pixelsize^2;
            base(i).(fn).area      = cellfun(@str2double, S.object_area(I))      * pixelsize^2;
            base(i).(fn).perimferet = cellfun(@str2double, S.object_feret(I))    * pixelsize;
            base(i).(fn).ESD       = 2 * sqrt((cellfun(@str2double, S.object_area(I)) * pixelsize^2) / pi);
        end

        % Conversion factor [m-3]
        base(i).(fn).conver = base(i).(fn).Fracnb ./ base(i).vol;

    end
end

close(h);

m         = numel(base);
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
