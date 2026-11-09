# FlowCam

The FlowCam is an imaging device for microplankton (Sieracki et al., 1998). Images are processed with VisualSpreadsheet and uploaded to EcoTaxa. In the pipeline, choose **Flowcam** in the instrument list.

Back to: [Getting started](../getting-started.md) · [Outputs](../outputs.md)

## 1. Required metadata

These fields are used in calculations: they must be filled in, in the expected units.

| EcoTaxa field | Content | Expected unit |
|---|---|---|
| `sample_id` | Sample identifier | |
| `sample_initial_col_vol_m3` | Volume of water collected (e.g. filtered by the net) | m³ |
| `sample_conc_vol_ml` (or `sample_volconc`) | Volume of the concentrated sample | mL |
| `acq_fluid_volume_imaged` | Volume of fluid imaged by the FlowCam | mL |
| `acq_raw_image_total` | Total number of raw images acquired | |
| `process_nb_images` | Number of raw images processed | |
| `process_pixel` | Pixel size | µm/pixel |
| `acq_id` | Size fraction identifier | |
| `object_lat`, `object_lon` | Sampling position | decimal degrees |
| `object_date`, `object_time` | Sampling date and time | |
| `object_major`, `object_minor`, `object_feret` | Major axis, minor axis, Feret diameter | pixels |
| `object_area`, `object_area_exc` | Area, area excluding holes | pixels |
| `object_annotation_category`, `object_annotation_hierarchy`, `object_annotation_status` | Taxonomic annotation and validation status | |

These fields are read and stored, but not used in calculations. They must be present in the export: `sample_ship`, `object_depth_min`, `object_depth_max`, `acq_min_esd`, `acq_max_esd`. The fields `object_lat_end` and `object_lon_end` are optional.

## 2. Volume conversion

The volume actually processed is first corrected for the images that were not processed:

$$V_{processed} = \frac{V_{imaged}}{N_{raw}} \times \min(N_{raw},\ N_{processed})$$

Then, for each fraction:

$$conver = \frac{V_{concentrated}}{V_{processed} \times V_{collected} \times f}$$

- $V_{imaged}$: volume of fluid imaged (`acq_fluid_volume_imaged`, converted from mL to m³);
- $N_{raw}$: total number of raw images acquired (`acq_raw_image_total`);
- $N_{processed}$: number of raw images processed (`process_nb_images`);
- $V_{concentrated}$: volume of the concentrated sample (`sample_conc_vol_ml`, converted from mL to m³);
- $V_{collected}$: volume of water collected (`sample_initial_col_vol_m3`, m³);
- $f$: fraction of the images annotated (1 by default; e.g. 0.3 if 30 % were annotated, see section 3).

Abundance (ind m⁻³) = number of objects × `conver`.

## 3. Instrument-specific steps: one or two fractions

After choosing the instrument, a window asks: **FlowCam images could be separated in two size fractions with different sub-sampling. Is that the case?**

**No (single fraction, most common case).** All the tsv files are in one folder. All objects are processed with $f = 1$.

**Yes (two fractions).** This option is used when only part of the images was annotated for one size fraction, to reduce the annotation effort. For example, in the Tara Pacific FlowCam datasets (Test datasets), objects > 45 µm were fully annotated (100 %), whereas a random subset of 30 % of the objects < 45 µm was annotated (Mériguet et al., 2025).

1. Export each fraction separately from EcoTaxa (one project or one export per fraction).
2. Put each fraction in its own folder. **The tsv files of a given sample must have the same name in both folders**: the pipeline pairs the two fractions of a sample by file name.
3. Select the first folder and enter its annotated percentage (e.g. 100), then the second folder and its percentage (e.g. 30).

The pipeline divides the conversion factor of each fraction by its annotated fraction $f$, so that abundances and biovolumes are representative of the whole sample.

## 4. Specificities

- All three biovolume methods are available (plain area, riddled area, ellipsoid).

## 5. Test datasets

**Single fraction**: folder `test_data/flowcam/single_fraction/`, 3 samples from Tonga expedition ([EcoTaxa project #2841](https://ecotaxa.obs-vlfr.fr/prj/2841)).

1. Run `EcoTaxa_Tool_Box`, choose **EcoTaxa export (tsv files) — process now**.
2. Select `test_data/flowcam/single_fraction/`, then **Flowcam**.
3. Answer **No** to the question about the two fractions, then follow [Getting started](../getting-started.md).

**Two fractions**: folder `test_data/flowcam/two_fractions/`, 3 samples from the Tara Pacific expedition : 100 % annotated ([EcoTaxa project #11341](https://ecotaxa.obs-vlfr.fr/prj/11341)) and 30 % annotated ([EcoTaxa project #11353](https://ecotaxa.obs-vlfr.fr/prj/11353)).

1. Run `EcoTaxa_Tool_Box`, choose **EcoTaxa export (tsv files) — process now**.
2. Select `test_data/flowcam/two_fractions/` (any folder: it is not used in this case), then **Flowcam**.
3. Answer **Yes**, select `fraction_100/` and enter 100, then `fraction_30/` and enter 30.
4. Follow [Getting started](../getting-started.md).