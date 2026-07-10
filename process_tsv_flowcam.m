% process_tsv_flowcam.m
% Assembles data from EcoTaxa TSV files (FlowCam) into the structured "base" database.

function [base, m, zoo_groups, Idlistshort] = process_tsv_flowcam

%% Optional two-fraction setup (100 % and subsampled fraction)

answer2 = questdlg( ...
    'FlowCam images could be separated in two size fractions with different sub-sampling. Is that the case?', ...
    'Fractionation', 'Yes', 'No', 'No');

folder1 = '';  folder2 = '';
answerfrac1 = 1;  answerfrac2 = 1;

switch answer2

    case 'Yes'
        f = msgbox('Please select the folder containing the *.tsv files of the first fraction (100 %)');
        uiwait(f);
        folder1 = uigetdir;
        answerfrac1 = str2double(cell2mat(inputdlg({'Fractionation (in %)'}, 'Fraction 1', 1, {'100'}))) / 100;

        f = msgbox('Please select the folder containing the *.tsv files of the second fraction');
        uiwait(f);
        folder2 = uigetdir;
        answerfrac2 = str2double(cell2mat(inputdlg({'Fractionation (in %)'}, 'Fraction 2', 1, {'30'}))) / 100;

        cd(folder1);

    case 'No'
        % stay in current directory

end

A         = dir('*.tsv');
filenames = {A.name};
m         = numel(filenames);

%% Initialise

base        = [];
Idlist      = {};
Idlistshort = {};

h = waitbar(0, 'Loading TSV files and building base...');

for i = 1:m

    waitbar(i/m, h, sprintf('Loading TSV files: %d of %d', i, m));

    switch answer2

        case 'Yes'
            cd(folder1);
            S = readtable(filenames{i}, 'FileType', 'text', 'ReadVariableNames', true);
            S.acq_id(:) = {'d1'};

            cd(folder2);
            S2 = readtable(filenames{i}, 'FileType', 'text', 'ReadVariableNames', true);
            S2.acq_id(:) = {'d2'};

            % Harmonise object_link type before concatenation
            if isnumeric(S.object_link),  S.object_link  = num2cell(S.object_link);  end
            if isnumeric(S2.object_link), S2.object_link = num2cell(S2.object_link); end
            S = [S; S2];

        case 'No'
            S = readtable(filenames{i}, 'FileType', 'text', 'ReadVariableNames', true);
    end

    %% Clean numeric columns stored as mixed cell arrays
    varnames = S.Properties.VariableNames;
    for j = 21:numel(varnames)
        col = S.(varnames{j});
        if iscellstr(col)
            emptyIdx = cellfun(@isempty, col);
            if any(emptyIdx)
                col(emptyIdx) = {'NaN'};
                S.(varnames{j}) = cellfun(@str2double, col);
            end
        end
    end

    %% Homogenise taxonomic annotation hierarchy
    S.object_annotation_hierarchy = strrep(S.object_annotation_hierarchy, '-', '_');
    S.object_annotation_hierarchy = strrep(S.object_annotation_hierarchy, '>', '_');

    sample = unique(S.sample_id);

    %% Build base — metadata
    base(i).sampleid  = sample;
    base(i).ship      = unique(S.sample_ship);
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

    base(i).zmax = unique(S.object_depth_max);
    base(i).zmin = unique(S.object_depth_min);

    %% Pixel size [mm/pixel]
    rawpix = unique(S.process_pixel);
    if isnumeric(rawpix)
        pixelsize = rawpix / 1000;
    else
        pixelsize = str2double(cell2mat(rawpix)) / 1000;
    end
    base(i).pixelsize = pixelsize;

    %% Volume fields — raw image count
    rawImgVal = unique(S.acq_raw_image_total);
    base(i).raw_image_total = toDouble(rawImgVal);

    nbImg = unique(S.process_nb_images);
    base(i).process_nb_images = toDouble(nbImg);

    %% Concentrated sample volume [m3]
    if ismember('sample_conc_vol_ml', S.Properties.VariableNames)
        base(i).volconc = toDouble(unique(S.sample_conc_vol_ml)) / 1e6;
    else
        base(i).volconc = toDouble(unique(S.sample_volconc)) / 1e6;
    end

    %% Initial collected volume [m3]
    if ismember('sample_initial_col_vol_m3', S.Properties.VariableNames)
        val = toDouble(unique(S.sample_initial_col_vol_m3));
    else
        val = toDouble(unique(S.sample_comment_or_volume)) / 1e6;
    end
    if isempty(val) || all(isnan(val)), val = NaN; end
    base(i).sample_initial_col_vol = val(1);

    %% Fluid volume imaged [m3]
    rawFluVol = unique(S.acq_fluid_volume_imaged);
    if isnumeric(rawFluVol)
        base(i).fluid_vol_imaged = rawFluVol / 1e6;
    else
        tmp = char(rawFluVol);
        idx = strfind(tmp, '_');
        base(i).fluid_vol_imaged = str2double(tmp(1:idx(1)-1)) / 1e6;
    end

    % Final processed volume [m3]
    base(i).volume_imaged_processed = ...
        (base(i).fluid_vol_imaged / base(i).raw_image_total) * ...
        min(base(i).raw_image_total, base(i).process_nb_images);

    %% Size fractions
    base(i).fracids = unique(S.acq_id);
    nfrac           = numel(base(i).fracids);

    for fracnb = 1:nfrac

        fn = ['d' num2str(fracnb)];
        I  = strcmp(S.acq_id, base(i).fracids{fracnb});

        base(i).(fn).fracmin  = unique(S.acq_min_esd);
        base(i).(fn).fracsup  = unique(S.acq_max_esd);
        base(i).(fn).anotStatus = S.object_annotation_status(I);

        % Taxonomic annotation
        base(i).(fn).object_annotation_hierarchy = S.object_annotation_hierarchy(I);

        % Update global taxonomic lists
        [Idlist, ia, ~] = unique([Idlist; S.object_annotation_hierarchy(I)]);
        Idlistshort     = [Idlistshort; S.object_annotation_category(I)];
        Idlistshort     = Idlistshort(ia);

        % Percent validated
        isVal = ~cellfun(@isempty, strfind(S.object_annotation_status(I), 'validated'));
        base(i).(fn).percentValidated = 100 * sum(isVal) / numel(isVal);

        % Morphological measurements
        if isnumeric(S.object_major)
            base(i).(fn).major       = S.object_major(I)    * pixelsize;
            base(i).(fn).minor       = S.object_minor(I)    * pixelsize;
            base(i).(fn).area_exc    = S.object_area_exc(I) * pixelsize^2;
            base(i).(fn).area        = S.object_area(I)     * pixelsize^2;
            base(i).(fn).area_origin = S.object_area(I);
            base(i).(fn).perimferet  = S.object_feret(I)    * pixelsize;
            base(i).(fn).ESD         = 2 * sqrt((S.object_area(I) * pixelsize^2) / pi);
        else
            base(i).(fn).major       = cellfun(@str2double, S.object_major(I))    * pixelsize;
            base(i).(fn).minor       = cellfun(@str2double, S.object_minor(I))    * pixelsize;
            base(i).(fn).area_exc    = cellfun(@str2double, S.object_area_exc(I)) * pixelsize^2;
            base(i).(fn).area        = cellfun(@str2double, S.object_area(I))     * pixelsize^2;
            base(i).(fn).area_origin = cellfun(@str2double, S.object_area(I));
            base(i).(fn).perimferet  = cellfun(@str2double, S.object_feret(I))    * pixelsize;
            base(i).(fn).ESD         = 2 * sqrt((cellfun(@str2double, S.object_area(I)) * pixelsize^2) / pi);
        end

        % Conversion factor [m-3]
        switch answer2
            case 'Yes'
                falsemotoda = answerfrac1 * (fracnb == 1) + answerfrac2 * (fracnb == 2);
                base(i).(fn).conver = base(i).volconc ./ ...
                    (base(i).volume_imaged_processed .* base(i).sample_initial_col_vol .* falsemotoda);
            case 'No'
                base(i).(fn).conver = base(i).volconc ./ ...
                    (base(i).volume_imaged_processed .* base(i).sample_initial_col_vol);
        end

    end
end

close(h);

m          = numel(base);
zoo_groups = Idlist;

end

%% Helper — convert numeric or cell-string scalar to double
function v = toDouble(x)
    if isnumeric(x)
        v = x;
    else
        v = str2double(cell2mat(x));
    end
end
