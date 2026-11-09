# EcoTaxa Tools Box Pipeline (MATLAB)

The EcoTaxa Tools Box is an open-source pipeline designed to process and explore taxonomically annotated images hosted on [EcoTaxa](https://ecotaxa.obs-vlfr.fr/). 

> An equivalent version is available in R: at [HZaccomer/EcotaxaTools](https://github.com/HZaccomer/EcotaxaTools).

It reads the tsv files exported from EcoTaxa for one instrument, merges them into a single database and computes key ecological descriptors for each sample. Each instrument has specific variable names and units: the pipeline accounts for them, so that the descriptors are computed in the same way and expressed in the same units whatever the instrument (IFCB, FlowCam, PlanktoScope, ZooScan). For each sample, it computes:

- abundance (ind m⁻³) and biovolume (mm³ m⁻³), per taxonomic annotation and per grouping level (living vs non-living, plankton groups, trophic groups);
- normalized biovolume size spectra (NBSS), with slope and intercept;
- summary figures per station and for the whole dataset (maps, compositions, size spectra, trophic pyramids).

## How it works

Run a single script: **`EcoTaxa_Tool_Box.m`**. It guides you through two successive parts:

1. **Processing**: the tsv files exported from EcoTaxa are merged into a single database, and abundances, biovolumes and NBSS are computed. This part is run automatically by `EcoTaxa_Tool_Box.m` (it calls `ecotaxa_tools.m`; you do not need to open it). The result is saved as `base_<instrument>_<project>.mat`.
2. **Analysis**: the `.mat` base is loaded and the figures are produced and saved as PDF, sorted into folders.

At launch, you choose where to start:

| Your data | Choice at launch | What happens |
|---|---|---|
| EcoTaxa tsv exports, not yet processed | **EcoTaxa export (tsv files) — process now** | Processing, then you can stop or continue to the analysis |
| A `.mat` base already created | **Matlab base.mat — already processed** | Analysis only, without reprocessing |

> An optional module for multivariate analyses (PCA, PCoA, hierarchical clustering) can then be run on the `.mat` base: `multivariates_analysis_following.m`.

## Supported instruments

| Instrument | Size range | Instrument-specific notes |
|---|---|---|
| [ZooScan](docs/instruments/zooscan.md) | > 200 µm | Optional several size fractions per sample (d1, d2), sub-sampling via `acq_sub_part` |
| [FlowCam](docs/instruments/flowcam.md) | 20–200 µm | Optional two fractions with different sub-sampling rates |
| [IFCB](docs/instruments/ifcb.md) | 10–150 µm | No `area_exc` field: no riddled-area biovolume |
| [PlanktoScope](docs/instruments/planktoscope.md) | 20–200 µm | Concentration and dilution factors required |

Each instrument page lists the mandatory EcoTaxa metadata, their expected units and the volume conversion formula.

## Requirements

- MATLAB R2025b (tested). Earlier versions: not guaranteed.
- No additional MATLAB toolbox required.
- Third-party libraries, included in `lib/`: [Fathom](https://www.usf.edu/marine-science/research/matlab-resources/fathom-toolbox-for-matlab.aspx), [M_Map](https://www.eoas.ubc.ca/~rich/map.html), `settingsdlg`.

## Installation

Clone or download this repository:

```
git clone https://github.com/zoemgt/EcoTaxa-Tools-Box-Pipeline.git
```

In MATLAB, set the repository as the current folder and add it to the path:

```matlab
addpath(genpath(pwd))
```

## Quick start (with the test datasets)

A small test dataset (3 samples) is provided for each instrument in `test_data/`:

| Instrument | Folder | Source |
|---|---|---|
| ZooScan | `test_data/zooscan/` | Tara Pacific expedition (Mériguet et al., 2025), [EcoTaxa project #11292](https://ecotaxa.obs-vlfr.fr/prj/11292) |
| FlowCam (single fraction) | `test_data/flowcam/single_fraction/` | [to be completed] |
| FlowCam (two fractions) | `test_data/flowcam/two_fractions/` | Tara Pacific expedition (Mériguet et al., 2024, 2025): 100 % annotated ([EcoTaxa project #11341](https://ecotaxa.obs-vlfr.fr/prj/11341)) and 30 % annotated ([EcoTaxa project #11353](https://ecotaxa.obs-vlfr.fr/prj/11353)) |
| IFCB | `test_data/ifcb/` | [to be completed] |
| PlanktoScope | `test_data/planktoscope/` | Ada-II 2023 expedition (Mériguet et al., 2026, under review), [EcoTaxa project #10056](https://ecotaxa.obs-vlfr.fr/prj/10056) |

1. In MATLAB, run `EcoTaxa_Tool_Box`.
2. Choose **EcoTaxa export (tsv files) — process now**.
3. Select the test folder of your instrument (e.g. `test_data/planktoscope/`), then the matching instrument in the list.
4. Answer the dialogs (taxonomic grouping, file name), then save the base.
5. Choose **Continue to ecological analysis**: the figures are written to the folders `global raw analysis/`, `NBSS analysis/` and to one sheet per station.

Each dialog is explained step by step in [Getting started](docs/getting-started.md), and instrument-specific steps in the instrument pages.

## Before using your own data

- Export from EcoTaxa with **Export → General export**, and set *Separate (in ZIP sub-directories) output by* to **Sample**.
- Check that all mandatory metadata are filled in and in the expected units (see the instrument pages). Missing or badly formatted metadata stop the pipeline or produce wrong values.

## Documentation

- [Getting started](docs/getting-started.md): step-by-step walkthrough of both parts
- [Outputs](docs/outputs.md): structure of the `.mat` base, variables and units
- Instrument pages: [ZooScan](docs/instruments/zooscan.md) · [FlowCam](docs/instruments/flowcam.md) · [IFCB](docs/instruments/ifcb.md) · [PlanktoScope](docs/instruments/planktoscope.md)

## Citation

If you use this pipeline, please cite:

Mériguet, Z., Perruchon, A., Zaccomer, H., et al. (2026) Studying plankton communities from annotated quantitative images on EcoTaxa: an automated plankton community processing and analysis pipeline. *Journal of Plankton Research* (under review).

Software version: [DOI Zenodo – to be added]

## License and contact

MIT License. See [LICENSE](LICENSE).

Zoé Mériguet (zoemeriguet@gmail.com), Fabien Lombard (fabien.lombard@imev-mer.fr), Laboratoire d'Océanographie de Villefranche-sur-Mer.
 
