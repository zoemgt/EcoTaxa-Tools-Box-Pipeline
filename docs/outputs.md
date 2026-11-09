# Outputs

Part 1 of the pipeline saves a MATLAB file `base_<instrument>_<project>.mat`. It contains a single variable, `base`: a structure array with one element per sample (`base(1)`, `base(2)`…). Part 2 and the multivariate module read this file.

```matlab
load('base_<instrument>_<project>.mat')   % replace with your file name; loads the structure 'base'
numel(base)                               % number of samples
base(1)                                   % fields of the first sample
```

Each sample contains four levels of information:

| Level | Content |
|---|---|
| Sample metadata | Identifier, date, position, volumes, pixel size… |
| `d1`, `d2`… | Results per size fraction (as defined in EcoTaxa `acq_id`) |
| `tot` | Results for the whole sample (sum of the fractions) |
| `regroupped` | Results per taxonomic annotation and per grouping level, NBSS and NBSS slope |

## 1. Size classes

All descriptors are computed on a common grid of **213 biovolume size classes**, identical for all instruments. The classes follow a geometric progression: the upper limit of each class is 2^0.25 times its lower limit (Platt & Denman, 1978). The grid covers biovolumes from 10⁻¹² to about 10⁴ mm³, i.e. equivalent spherical diameters (ESD) from about 0.1 µm to 27 mm. Each instrument only fills the classes corresponding to its own size range.

| Variable | Content | Unit |
|---|---|---|
| `X` | Midpoint of each size class (natural logarithm) | ln(mm³) |
| `X1` | Width of each size class, used to de-normalize NBSS | mm³ |
| `tot.ESDvector` | Midpoint of each size class expressed as ESD | µm |

Note that `X` is stored as a natural logarithm. To use it with log10-transformed values (e.g. to fit a NBSS slope), convert it first: `X / log(10)`.

## 2. Sample metadata

Common to all instruments:

| Variable | Content |
|---|---|
| `sampleid` | EcoTaxa sample identifier |
| `date`, `time` | Sampling date and time |
| `latitude`, `longitude` | Sampling position (start) |
| `pixelsize` | Pixel size (mm/pixel) |
| `fracids` | Identifiers of the size fractions (`acq_id`) |

Other metadata (volumes, end position, depth, net type, ship…) depend on the instrument: see the instrument pages.

## 3. Results per fraction (`d1`, `d2`…) and per sample (`tot`)

For each fraction, `base(i).d1`, `base(i).d2`… contain the morphometric measurements of each object (`major`, `minor` in mm; `area`, `area_exc` in mm²; `ESD`), its taxonomic annotation, and the following descriptors. `base(i).tot` contains the same descriptors summed over all fractions.

| Variable | Content | Unit |
|---|---|---|
| `conver` | Conversion factor from imaged objects to concentration (fractions only) | m⁻³ |
| `Ab` | Abundance per taxonomic annotation | ind m⁻³ |
| `Abtot` | Number of imaged objects per taxonomic annotation (fractions only) | objects |
| `Bv` | Biovolume per taxonomic annotation, ellipsoidal method | mm³ m⁻³ |
| `Bvtot` | Imaged biovolume per taxonomic annotation, not converted (fractions only) | mm³ |
| `Yab` | Abundance size spectrum (normalized by class width) | ind m⁻³ mm⁻³ |
| `Ybv_Plain_Area_BV_spectra` | NBSS, biovolume from object area | mm³ m⁻³ mm⁻³ |
| `Ybv_Riddled_Area_BV_spectra` | NBSS, biovolume from excluded area (not available for IFCB) | mm³ m⁻³ mm⁻³ |
| `Ybv_Ellipsoid_BV_spectra` | NBSS, ellipsoidal biovolume (from major and minor axes) | mm³ m⁻³ mm⁻³ |
| `Zoo_groups` | Taxonomic annotations, in the order of the columns | |

Spectra are matrices of 213 rows (size classes) × one column per taxonomic annotation.

The three biovolume methods are described in Vandromme et al. (2012). The ellipsoidal method is recommended: it is used for `Bv`, for the NBSS of section 4 and for all figures.

## 4. Grouped results (`regroupped`)

`base(i).regroupped` contains the same descriptors (`Ab`, `Bv`, `Yab`, `Ybv_*`), but the columns are arranged in successive blocks:

| Block | Content | Position |
|---|---|---|
| Taxonomic annotations | One column per taxon | `originalplace` (first and last column) |
| All | Sum of all objects | `allplace` |
| Living / non-living | Living, then non-living | `livingplace`, `notlivingplace` |
| Plankton groups | One column per plankton group | `planktongroupplace` (first and last column) |
| Trophic groups | One column per trophic group code | `trophicplace` (first and last column) |

The label of each column is given in `regroupped.Zoo_groups`. Trophic groups are labelled with their code (see [Getting started](getting-started.md#4-taxonomic-grouping)).

`regroupped` also contains the NBSS used in part 2, and its fit:

| Variable | Content |
|---|---|
| `NBSSliving`, `NBSSpft`, `NBSStrophic` | NBSS (ellipsoidal biovolume) of living organisms, per plankton group, per trophic group |
| `NBSSliving_trunked`, `NBSSpft_trunked`, `NBSStrophic_trunked` | Same, truncated: classes below the spectrum maximum, empty classes and isolated tails are removed |
| `NBSSliving_slope`, `NBSSliving_intercept` | Slope and intercept of the linear fit log10(NBSS) = a + b · log10(biovolume), on the truncated NBSS of living organisms. The intercept is log10(NBSS) at a biovolume of 1 mm³. NaN when fewer than 5 size classes are available |

## 5. Interpreting the outputs

- Abundances and biovolumes depend on the quality of the annotation and of the metadata: check volumes and units before any interpretation.
- Size classes at the edges of the instrument's size range are under-sampled: the truncation of the NBSS removes the smallest classes for this reason.
- The figures produced by part 2 are exploratory. Choose the descriptors and methods adapted to your scientific question before any publication.