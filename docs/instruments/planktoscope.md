# PlanktoScope

The PlanktoScope is a low-cost, open-source imaging platform for microplankton (Pollina et al., 2022). Images are segmented by the PlanktoScope software and uploaded to EcoTaxa. In the pipeline, choose **PlanktoScope** in the instrument list.

Back to: [Getting started](../getting-started.md) · [Outputs](../outputs.md)

## 1. Required metadata

These fields are used in calculations: they must be filled in, in the expected units.

| EcoTaxa field | Content | Expected unit |
|---|---|---|
| `sample_id` | Sample identifier | |
| `sample_total_volume` | Volume of water collected (e.g. filtered by the net) | L |
| `sample_concentrated_sample_volume` | Volume of the concentrated sample | mL |
| `sample_dilution_factor` or `sample_concentration_factor` | Dilution or concentration factor (one of the two; 1 if both are empty, see section 2) | no unit |
| `acq_imaged_volume` | Volume imaged by the PlanktoScope | mL |
| `process_pixel` | Pixel size | µm/pixel |
| `acq_id` | Acquisition identifier | |
| `object_lat`, `object_lon` | Sampling position | decimal degrees |
| `object_date`, `object_time` | Sampling date and time | |
| `object_major`, `object_minor` | Major and minor axes | pixels |
| `object_area`, `object_area_exc` | Area, area excluding holes | pixels |
| `object_annotation_category`, `object_annotation_hierarchy`, `object_annotation_status` | Taxonomic annotation and validation status | |

Optional: `object_lat_end`, `object_lon_end`.

These fields are read and stored, but not used in calculations. They must be present in the export: `sample_ship`, `sample_project`, `sample_sampling_gear`, `object_depth_min`, `object_depth_max`, `acq_nb_frame`, `acq_celltype`, `acq_camera_resolution`.

## 2. Volume conversion

$$conver = \frac{V_{concentrated} \times F}{V_{imaged} \times V_{total}}$$

- $V_{concentrated}$: volume of the concentrated sample (`sample_concentrated_sample_volume`, converted from mL to m³);
- $V_{imaged}$: volume imaged (`acq_imaged_volume`, converted from mL to m³);
- $V_{total}$: volume of water collected (`sample_total_volume`, converted from L to m³);
- $F$: dilution or concentration correction, chosen in this order:
  1. if `sample_dilution_factor` is filled in: $F$ = `sample_dilution_factor` (e.g. 2 if the sample was diluted twice before imaging);
  2. otherwise, if `sample_concentration_factor` is filled in: $F$ = 1 / `sample_concentration_factor`;
  3. otherwise: $F$ = 1 (no correction).

Fill in only one of the two fields. The value used for each sample is stored in the base (`volume_factor`).

Abundance (ind m⁻³) = number of objects × `conver`.

## 3. Instrument-specific steps

No additional window.

## 4. Specificities

- All three biovolume methods are available (plain area, riddled area, ellipsoid).
- Small, round objects close to the resolution limit of the PlanktoScope are difficult to identify: they are often annotated as unidentified.

## 5. Test dataset

Folder `test_data/planktoscope/`: 3 samples from the Ada-II 2023 expedition (Mériguet et al., 2026, under review), [EcoTaxa project #10056](https://ecotaxa.obs-vlfr.fr/prj/10056).

1. Run `EcoTaxa_Tool_Box`, choose **EcoTaxa export (tsv files) — process now**.
2. Select `test_data/planktoscope/`, then **PlanktoScope**.
3. Follow the steps described in [Getting started](../getting-started.md).