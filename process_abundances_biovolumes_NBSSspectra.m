% process_abundances_biovolumes_NBSSspectra.m
% Computes abundance and biovolume (3 methods) for all annotated groups,
% plus normalised biovolume size spectra (NBSS).
%
% Units:
%   abundance  : ind/m3
%   biovolume  : mm3/m3
%   spectra Y  : m3  (log-transform before plotting)

function [base_spectres, SStot] = process_abundances_biovolumes_NBSSspectra( ...
    subset, base_spectres, smin, smax, k, uu, zoo_groups, i, fracnb)

isTot = strcmp(fracnb, 'tot');
SStot = [];

% Clean annotation hierarchy in group list
zoo_groups = strrep(zoo_groups, '-', '_');

m = numel(zoo_groups);

%% Size-class bin vector
xvec = smin;
while xvec(end) < smax
    xvec(end+1, 1) = smin * k^(numel(xvec)); %#ok<AGROW>
end
X1   = diff(xvec);
xlog = log(xvec);
X    = (xlog(1:end-1) + xlog(2:end)) / 2;  % bin centres (log scale)
nb_int = numel(X);

%% Extract morphological data
mj       = subset.major;    % mm
mn       = subset.minor;    % mm
aarea    = subset.area;     % mm2
aareaexc = subset.area_exc; % mm2
conver   = subset.conver;   % m-3

% Store conversion factor
fn = fieldName(isTot, fracnb);
base_spectres(i).(fn).conver = conver;

%% Build taxonomic identity matrix once (shared across methods)
pred = subset.object_annotation_hierarchy;
ID   = false(numel(pred), m);
for j = 1:m
    ID(:,j) = strcmp(zoo_groups{j}, pred);
end
ID = double(ID);

%% Loop over 3 biovolume methods
method_list = {'Plain_Area_BV_spectra', 'Riddled_Area_BV_spectra', 'Ellipsoid_BV_spectra'};

for meth = 1:3

    % Biovolume calculation [mm3]
    switch meth
        case 3  % Ellipsoid (recommended)
            sss = (4/3) * pi .* (mj/2) .* (mn/2) .* (mn/2);
        case 2  % Riddled area
            esd = 2 * sqrt(aareaexc / pi);
            sss = (4/3) * pi .* (esd/2).^3;
        case 1  % Plain area
            esd = 2 * sqrt(aarea / pi);
            sss = (4/3) * pi .* (esd/2).^3;
    end

    vol = double(sss);
    sss = double(log(sss));

    %% Per-group abundance, biovolume and spectra
    Ab    = zeros(m, 1);
    Abtot = zeros(m, 1);
    Bv    = zeros(m, 1);
    Bvtot = zeros(m, 1);
    Yab   = zeros(nb_int, m);
    Ybv   = zeros(nb_int, m);

    for j = 1:m

        idx_j = logical(ID(:,j));

        % Abundance
        Ab(j)    = sum(idx_j) .* conver;
        Abtot(j) = sum(idx_j);

        % Biovolume (Ellipsoid method only, stored once)
        if meth == 3
            Bv(j)    = sum(vol(idx_j)) .* conver;
            Bvtot(j) = sum(vol(idx_j));
        end

        % Size spectra
        if any(idx_j)
            [~, Yab_j, ~, ~] = f_ZOO_spectrum(sss(idx_j), vol(idx_j), smin, smax, k, 1);
            [~, Ybv_j, ~, ~] = f_ZOO_spectrum(sss(idx_j), vol(idx_j), smin, smax, k, 2);
            Yab(:,j) = Yab_j .* conver;
            Ybv(:,j) = Ybv_j .* conver;
        end
        % zeros already assigned for absent groups

    end

    %% Store in base_spectres
    base_spectres(i).(fn).Ab    = Ab;
    base_spectres(i).(fn).Abtot = Abtot;

    if meth == 3
        base_spectres(i).(fn).Bv    = Bv;
        base_spectres(i).(fn).Bvtot = Bvtot;
    end

    base_spectres(i).(fn).Yab = Yab;
    base_spectres(i).(fn).(['Ybv_' method_list{meth}]) = Ybv;

end

base_spectres(i).(fn).X          = X;
base_spectres(i).(fn).X1         = X1;
base_spectres(i).(fn).Zoo_groups = zoo_groups;

end

%% Helper — return struct field name for tot vs d1/d2
function fn = fieldName(isTot, fracnb)
    if isTot
        fn = 'tot';
    else
        fn = ['d' num2str(fracnb)];
    end
end
