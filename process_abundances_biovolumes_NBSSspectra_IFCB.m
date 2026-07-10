% process_abundances_biovolumes_NBSSspectra_IFCB.m
% Computes abundance and biovolume (3 methods) for all annotated IFCB groups,
% plus normalised biovolume size spectra (NBSS).
%
% Units:
%   abundance  : ind/m3
%   biovolume  : mm3/m3
%   spectra Y  : m3  (log-transform before plotting)
%
% IFCB-specific method 3 uses the pre-computed summed biovolume field.

function [base_spectres, SStot] = process_abundances_biovolumes_NBSSspectra_IFCB( ...
    subset, base_spectres, smin, smax, k, uu, zoo_groups, i, fracnb)

isTot = strcmp(fracnb, 'tot');
SStot = [];

% Clean annotation hierarchy in group list
zoo_groups = strrep(zoo_groups, '-', '_');

m = numel(zoo_groups);

%% Size-class bin vector
xvec = smin;
while xvec(end) < smax
    xvec(end+1,1) = smin * k^(numel(xvec)); %#ok<AGROW>
end
X1   = diff(xvec);
xlog = log(xvec);
X    = (xlog(1:end-1) + xlog(2:end)) / 2;
nb_int = numel(X);

%% Extract morphological data
mj              = subset.major;           % mm
mn              = subset.minor;           % mm
aarea           = subset.area;            % mm2
summedbiovolume = subset.summedbiovolume; % mm3
conver          = subset.conver;          % m-3

fn = fieldName(isTot, fracnb);
base_spectres(i).(fn).conver = conver;

%% Taxonomic identity matrix
pred = subset.object_annotation_hierarchy;
ID   = false(numel(pred), m);
for j = 1:m
    ID(:,j) = strcmp(zoo_groups{j}, pred);
end
ID = double(ID);

%% 3 biovolume methods
method_list = {'Plain_Area_BV_spectra', 'Ellipsoid_BV_spectra', 'Summed_plain_Area_BV_spectra'};

for meth = 1:3

    switch meth
        case 3  % Summed biovolume from IFCB
            sss = summedbiovolume;
        case 2  % Ellipsoid
            sss = (4/3) * pi .* (mj/2) .* (mn/2) .* (mn/2);
        case 1  % Plain area
            esd = 2 * sqrt(aarea / pi);
            sss = (4/3) * pi .* (esd/2).^3;
    end

    vol = double(sss);
    vol(isnan(vol)) = 0;
    sss = double(log(sss));

    Ab    = zeros(m, 1);
    Abtot = zeros(m, 1);
    Bv    = zeros(m, 1);
    Bvtot = zeros(m, 1);
    Yab   = zeros(nb_int, m);
    Ybv   = zeros(nb_int, m);

    for j = 1:m

        idx_j = logical(ID(:,j));

        Ab(j)    = sum(idx_j) .* conver;
        Abtot(j) = sum(idx_j);

        if meth == 3
            Bv(j)    = sum(vol(idx_j)) .* conver;
            Bvtot(j) = sum(vol(idx_j));
        end

        if any(idx_j)
            [~, Yab_j, ~, ~] = f_ZOO_spectrum(sss(idx_j), vol(idx_j), smin, smax, k, 1);
            [~, Ybv_j, ~, ~] = f_ZOO_spectrum(sss(idx_j), vol(idx_j), smin, smax, k, 2);
            Yab(:,j) = Yab_j .* conver;
            Ybv(:,j) = Ybv_j .* conver;
        end

    end

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

function fn = fieldName(isTot, fracnb)
    if isTot
        fn = 'tot';
    else
        fn = ['d' num2str(fracnb)];
    end
end
