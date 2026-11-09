# ZooScan

The ZooScan is an automated imaging system for mesozooplankton (> 200 µm) (Gorsky et al., 2010). Images are segmented with ZooProcess and uploaded to EcoTaxa. In the pipeline, choose **Zooscan** in the instrument list.

Back to: [Getting started](../getting-started.md) · [Outputs](../outputs.md)

## 1. Required metadata

These fields are used in calculations: they must be filled in, in the expected units.

| EcoTaxa field | Content | Expected unit |
|---|---|---|
| `sample_id` | Sample identifier | |
| `sample_tot_vol` | Volume of water filtered by the net | m³ |
| `acq_id` | Size fraction identifier (e.g. d1, d2) | |
| `acq_sub_part` | Sub-sampling factor of the fraction (e.g. 4 if 1/4 of the sample was scanned) | |
| `process_particle_pixel_size_mm` | Pixel size | mm/pixel |
| `object_lat`, `object_lon` | Sampling position | decimal degrees |
| `object_date`, `object_time` | Sampling date and time | |
| `object_major`, `object_minor`, `object_feret` | Major axis, minor axis, Feret diameter | pixels |
| `object_area`, `object_area_exc` | Area, area excluding holes | pixels |
| `object_annotation_category`, `object_annotation_hierarchy`, `object_annotation_status` | Taxonomic annotation and validation status | |

These fields are read and stored, but not used in calculations. They must be present in the export and can be empty: `sample_scan_operator`, `sample_ship`, `sample_program`, `sample_stationid`, `sample_bottomdepth`, `sample_ctdrosettefilename`, `sample_other_ref`, `sample_barcode`, `sample_net_type`, `sample_tow_nb`, `sample_tow_type`, `sample_net_mesh`, `sample_net_surf`, `sample_zmin`, `sample_zmax`, `sample_comment`, `acq_min_mesh`, `acq_max_mesh`, `process_img_resolution`.

## 2. Volume conversion

For each fraction, the conversion factor from imaged objects to concentration is:

$$conver = \frac{acq\_sub\_part}{sample\_tot\_vol}$$

- `acq_sub_part`: sub-sampling factor of the fraction (number of parts the sample was split into);
- `sample_tot_vol`: volume filtered by the net (m³).

Abundance (ind m⁻³) = number of objects × `conver`.

## 3. Instrument-specific steps

No additional window. Size fractions are read from `acq_id`: each fraction is processed separately (`d1`, `d2`…), then summed in `tot` (see [Outputs](../outputs.md)). Check that the sub-sampling factor `acq_sub_part` is correct for each fraction.

## 4. Specificities

- All three biovolume methods are available (plain area, riddled area, ellipsoid).
- Touching objects must be separated in ZooProcess before upload: they bias abundance and size estimates (Vandromme et al., 2012).

## 5. Test dataset

Folder `test_data/zooscan/`: 3 samples from the Tara Pacific expedition : [EcoTaxa project #11292](https://ecotaxa.obs-vlfr.fr/prj/11292).

1. Run `EcoTaxa_Tool_Box`, choose **EcoTaxa export (tsv files) — process now**.
2. Select `test_data/zooscan/`, then **Zooscan**.
3. Follow the steps described in [Getting started](../getting-started.md).