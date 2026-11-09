# IFCB

The Imaging FlowCytobot (IFCB) is a submersible imaging flow cytometer for nano- and microplankton (Sosik & Olson, 2007). Images are processed with the IFCB software and uploaded to EcoTaxa. In the pipeline, choose **IFCB** in the instrument list.

Back to: [Getting started](../getting-started.md) · [Outputs](../outputs.md)

## 1. Required metadata

These fields are used in calculations: they must be filled in, in the expected units.

| EcoTaxa field | Content | Expected unit |
|---|---|---|
| `sample_id` | Sample identifier | |
| `acq_volume_sampled` | Volume of water analysed | mL |
| `acq_resolution_pixel_per_micron` | Image resolution | pixels/µm |
| `acq_id` | Acquisition identifier | |
| `object_lat`, `object_lon` | Sampling position | decimal degrees |
| `object_date`, `object_time` | Sampling date and time | |
| `object_major_axis_length`, `object_minor_axis_length` | Major and minor axes | pixels |
| `object_surface_area` | Area | pixels |
| `object_summed_biovolume` | Biovolume computed by the IFCB software | pixels³ |
| `object_summed_surface_area` | Summed surface area computed by the IFCB software | pixels |
| `object_annotation_category`, `object_annotation_hierarchy`, `object_annotation_status` | Taxonomic annotation and validation status | |

The fields `object_lat_end` and `object_lon_end` are optional.

## 2. Volume conversion

The IFCB analyses a known volume of water directly, without concentration:

$$conver = \frac{1}{V_{sampled}}$$

- $V_{sampled}$: volume of water analysed (`acq_volume_sampled`, converted from mL to m³).

Abundance (ind m⁻³) = number of objects × `conver`.

## 3. Instrument-specific steps

No additional window. The IFCB is processed by dedicated scripts (`process_spectra_and_save_IFCB.m`, `process_abundances_biovolumes_NBSSspectra_IFCB.m`), called automatically.

## 4. Specificities

- The IFCB export has no `area_exc` field: the riddled-area biovolume is not available.
- The three biovolume spectra available are: plain area (`Ybv_Plain_Area_BV_spectra`), ellipsoid (`Ybv_Ellipsoid_BV_spectra`) and the summed biovolume computed by the IFCB software (`Ybv_Summed_plain_Area_BV_spectra`).
- The biovolume per taxon `Bv` uses the summed biovolume computed by the IFCB software. The NBSS and figures use the ellipsoidal method, as for the other instruments.

## 5. Test dataset

Folder `test_data/ifcb/` ([to be completed]).

1. Run `EcoTaxa_Tool_Box`, choose **EcoTaxa export (tsv files) — process now**.
2. Select `test_data/ifcb/`, then **IFCB**.
3. Follow the steps described in [Getting started](../getting-started.md).