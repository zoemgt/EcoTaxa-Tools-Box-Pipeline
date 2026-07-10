% f_regroup_all.m
% Groups plankton abundance, biovolume, and size spectra by:
%   - original taxonomic annotations
%   - living / non-living
%   - plankton functional type (PFT) groups
%   - trophic groups
%
% Works for all instruments. The "Riddled Area" biovolume spectrum
% (Ybv_Riddled_Area_BV_spectra) is only present for instruments that export
% an "area_exc" field (ZooScan, FlowCam, PlanktoScope) — absent for IFCB,
% which is detected automatically via isfield() and skipped when missing.
%
% Input:
%   table_groupage : cell array with columns
%       1 = original taxon name
%       2 = living / non-living label
%       3 = PFT group name
%       4 = trophic level (numeric)
%   base           : single-sample struct (base(i).tot)
%
% Output:
%   base_regroup   : struct with grouped Ab, Bv, spectra, and group labels

function base_regroup = f_regroup_all(table_groupage, base)

tab      = table_groupage;
tab_orig = base.Zoo_groups;

%% Match original taxa to the groupage table (replaces strmatch)
to_keep = [];
for i = 1:numel(tab_orig)
    a = find(strcmp(tab_orig{i}, tab(:,1)));
    if ~isempty(a)
        if numel(a) > 1
            warning('f_regroup_all: "%s" appears multiple times in the groupage table.', tab_orig{i});
        end
        to_keep = [to_keep; a(1)];
    end
end

tab = tab(to_keep, :);

g_orig          = tab(:,1);
g_living        = tab(:,2);
g_planktongroup = tab(:,3);
g_trophic       = tab(:,4);

%% Extract base variables
hasRiddled = isfield(base, 'Ybv_Riddled_Area_BV_spectra');

Ab                        = base.Ab;
Ybv_Plain_Area_BV_spectra = base.Ybv_Plain_Area_BV_spectra;
Bv                        = base.Bv;
ss                        = base.Ybv_Ellipsoid_BV_spectra;
Yab                       = base.Yab;

if hasRiddled
    Ybv_Riddled_Area_BV_spectra = base.Ybv_Riddled_Area_BV_spectra;
end

%% Total (all groups summed)
Ab_t                        = sum(Ab);
Ybv_Plain_Area_BV_spectra_t = sum(Ybv_Plain_Area_BV_spectra, 2);
Bv_t                        = sum(Bv);
sst                         = sum(ss, 2);
Yab_t                       = sum(Yab, 2);

if hasRiddled
    Ybv_Riddled_Area_BV_spectra_t = sum(Ybv_Riddled_Area_BV_spectra, 2);
end

%% Grouping helper — builds summed arrays for any grouping variable
    function [Ab_n, Ybv_Plain_n, Ybv_Riddled_n, Bv_n, ssn, Yab_n, new_groups] = ...
            group_by(g_var, is_numeric_group)

        if nargin < 2, is_numeric_group = false; end

        if is_numeric_group
            new_groups = unique(cell2mat(g_var));
        else
            new_groups = unique(g_var);
        end

        ng       = numel(new_groups);
        Ab_n     = zeros(1, ng);
        Bv_n     = zeros(1, ng);
        ybv_rows = size(ss, 1);
        Ybv_Plain_n = zeros(ybv_rows, ng);
        ssn         = zeros(ybv_rows, ng);
        Yab_n       = zeros(ybv_rows, ng);
        if hasRiddled
            Ybv_Riddled_n = zeros(ybv_rows, ng);
        else
            Ybv_Riddled_n = [];
        end

        for gi = 1:ng
            if is_numeric_group
                f_ng = (new_groups(gi) == cell2mat(g_var));
            else
                f_ng = strcmp(new_groups{gi}, g_var);
            end
            Ab_n(gi)          = sum(Ab(f_ng));
            Bv_n(gi)          = sum(Bv(f_ng));
            Ybv_Plain_n(:,gi) = sum(Ybv_Plain_Area_BV_spectra(:,f_ng), 2);
            ssn(:,gi)         = sum(ss(:,f_ng), 2);
            Yab_n(:,gi)       = sum(Yab(:,f_ng), 2);
            if hasRiddled
                Ybv_Riddled_n(:,gi) = sum(Ybv_Riddled_Area_BV_spectra(:,f_ng), 2);
            end
        end
    end

%% Group by living / non-living
[Ab_n1, Ybv_Plain_n1, Ybv_Riddled_n1, Bv_n1, ssn1, Yab_n1, living_groups] = ...
    group_by(g_living);

%% Group by plankton functional type
[Ab_n2, Ybv_Plain_n2, Ybv_Riddled_n2, Bv_n2, ssn2, Yab_n2, pft_groups] = ...
    group_by(g_planktongroup);

%% Group by trophic level
[Ab_n3, Ybv_Plain_n3, Ybv_Riddled_n3, Bv_n3, ssn3, Yab_n3, trophic_groups] = ...
    group_by(g_trophic, true);

%% Assemble regrouped arrays
base_regroup = struct();

base_regroup.Ab = [Ab', Ab_t, Ab_n1, Ab_n2, Ab_n3];

base_regroup.Ybv_Plain_Area_BV_spectra = ...
    [Ybv_Plain_Area_BV_spectra, Ybv_Plain_Area_BV_spectra_t, ...
     Ybv_Plain_n1, Ybv_Plain_n2, Ybv_Plain_n3];

if hasRiddled
    base_regroup.Ybv_Riddled_Area_BV_spectra = ...
        [Ybv_Riddled_Area_BV_spectra, Ybv_Riddled_Area_BV_spectra_t, ...
         Ybv_Riddled_n1, Ybv_Riddled_n2, Ybv_Riddled_n3];
end

base_regroup.Bv = [Bv', Bv_t, Bv_n1, Bv_n2, Bv_n3];

base_regroup.Ybv_Ellipsoid_BV_spectra = ...
    [ss, sst, ssn1, ssn2, ssn3];

base_regroup.Yab = [Yab, Yab_t, Yab_n1, Yab_n2, Yab_n3];

%% Index positions for each grouping level
n_orig    = numel(g_orig);
n_pft     = numel(pft_groups);
n_trophic = numel(trophic_groups);

base_regroup.Zoo_groups = [g_orig; {'all'}; living_groups(:); pft_groups(:); ...
                            cellstr(num2str(trophic_groups(:)))];

base_regroup.originalplace     = [1, n_orig];
base_regroup.allplace           = n_orig + 1;
base_regroup.livingplace        = n_orig + 2;
base_regroup.notlivingplace     = n_orig + 3;
base_regroup.planktongroupplace = [n_orig + 4, n_orig + 3 + n_pft];
base_regroup.trophicplace       = [n_orig + 3 + n_pft + 1, ...
                                   n_orig + 3 + n_pft + n_trophic];

end
