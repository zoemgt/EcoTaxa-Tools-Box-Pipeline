# Getting started

This page walks you through a complete run of the EcoTaxa Tools Box, from the EcoTaxa export to the figures. The pipeline is interactive: each step is a dialog window. The tables below list the windows in the order in which they appear, and what to answer.

Instrument-specific windows and metadata are described in the instrument pages: [ZooScan](instruments/zooscan.md) · [FlowCam](instruments/flowcam.md) · [IFCB](instruments/ifcb.md) · [PlanktoScope](instruments/planktoscope.md).

## 1. Before you start

**Export your data from EcoTaxa**

1. In your EcoTaxa project, go to **Export → General export**.
2. Keep the default settings, except *Separate (in ZIP sub-directories) output by*, which must be set to **Sample**.
3. Unzip the archive and gather all the tsv files (one per sample) in a single folder.

**Check your metadata**

All mandatory metadata must be filled in, in the expected units. The list depends on the instrument: see the instrument pages. Missing or badly formatted metadata stop the pipeline or produce wrong values.

**First run?** Use a test dataset from `test_data/` (see the [README](../README.md#quick-start-with-the-test-datasets)).

## 2. Launch the pipeline

In MATLAB, set the repository folder as the current folder, add it to the path and run the main script:

```matlab
addpath(genpath(pwd))
EcoTaxa_Tool_Box
```

The first window asks **What data should be loaded?**

| Choice | When to use it | Next step |
|---|---|---|
| EcoTaxa export (tsv files) — process now | First run on a dataset | Part 1, then part 2 |
| Matlab base.mat — already processed | The `.mat` base already exists | Part 2 directly |

## 3. Part 1 — Processing the tsv files

| Window | What to do |
|---|---|
| *Before you start* (reminder) | Check that the export is "by sample" and the metadata complete, then click OK. |
| Select the folder containing the tsv files | Select the folder with your tsv files (one file per sample). |
| *Which instrument was used?* | Choose ZooScan, FlowCam, IFCB or PlanktoScope. |
| Instrument-specific windows | E.g. FlowCam size fractions. See the instrument page. |
| *Trophic group legend* | Information only: codes used for trophic groups (see section 4). |
| *Temporary categories mapping* (only if t00X categories are present) | See section 4. |
| *New taxa functional mapping* (only if new taxa are found) | See section 4. |
| *Save base as* | Enter the instrument and project names. The base is saved as `base_<instrument>_<project>.mat`. |
| *Next step?* | **Stop after processing**, or **Continue to ecological analysis** (part 2). If you continue, select the `.mat` base you have just saved. |

During processing, the pipeline computes for each sample the abundance, the biovolume (three methods) and the NBSS, per taxonomic annotation and per grouping level. The content of the base is described in [Outputs](outputs.md).

## 4. Taxonomic grouping

Each taxonomic annotation is assigned to three grouping levels, defined in the file `trophic_affiliation_of_organisms.xlsx` (repository root):

- **living / non-living**;
- **plankton group** (e.g. Bacillariophyta, Copepoda, detritus);
- **trophic group**, coded as follows:

| Code | Trophic group |
|---|---|
| -1 | not feeding |
| 1 | phototrophs |
| 1.5 | mixotrophs |
| 2 | grazers |
| 2.5 | omnivorous |
| 3 | predators |
| 3.5 | unknown |

You can edit this file to adapt the affiliations to your dataset.

**New taxa.** If your dataset contains taxa that are not yet in the file, a window opens for each of them: choose its living status, plankton group and trophic group.

> **Warning.** These choices are written permanently into `trophic_affiliation_of_organisms.xlsx` and will be applied to all your future runs, for every dataset. Check each choice carefully before validating. A backup copy of the previous file is saved automatically in the repository folder (`trophic_affiliation_of_organisms_backup_<date>.xlsx`): to undo a mistake, restore it or edit the file directly in Excel.

**Temporary categories (t00X).** Objects that cannot be identified are usually annotated in broad categories such as "other<plankton", which may end up containing a large and heterogeneous part of the plankton. As an alternative, used at the Quantitative Imaging Platform of Villefranche-sur-Mer (PIQv), visually consistent groups of unidentified objects can be annotated in EcoTaxa as temporary morphotypes (`t001`, `t002`…). This keeps them isolated: they can be identified later, or analysed separately.

Using these categories is optional. If your project does not contain any, the pipeline skips this step. If it does, their meaning is specific to your project, so they are not stored in the main affiliation file: the pipeline offers three options.

| Choice | What happens |
|---|---|
| Yes — load them | Select an existing mapping file (`.xlsx`) for these categories. |
| No — create them | A window opens for each category: choose its living status, plankton group and trophic group. The mapping is saved in a separate `.xlsx` file, which you can reload for the next runs. |
| No — ignore them (not recommended) | These objects are removed from the analysis. |

## 5. Part 2 — Analysis

| Window | What to do |
|---|---|
| Select the saved `*.mat` base file | Select your base (only if you started from an existing base). |
| *Export all graphs as PDF?* | **Yes** to save the figures; **No** to display them only. |
| *Exclude some taxonomic annotation from analysis?* | See below. |
| *Exclude some plankton groups from analysis?* | Same choices, at the plankton group level. |
| *Exclude some trophic groups from analysis?* | Same choices, at the trophic group level. |
| *Data for bar charts and map per station* | Abundance or biovolume. |
| *Select one taxon/group for dedicated analysis?* | Optional: maps and bar charts for one taxon or group, or for each of them. |
| *Run optional multivariate analyses?* | Optional: see section 6. |

**Excluding groups.** At each level, you can:

- **Exclude manually**: select the groups to remove (Ctrl/Cmd for multiple selection), for example detritus or other non-living objects to focus on living organisms. The list is saved as `<level>_excluded.xlsx`;
- **Load exclusion list**: reload a list saved during a previous run;
- **No exclusion**.

The NBSS analysis (spectra, truncated spectra, map of slopes) runs automatically.

## 6. Optional multivariate analyses

At the end of part 2, the pipeline offers multivariate analyses. This module is an example tool for exploring the data: normalization, ordination and clustering methods should be selected according to your dataset and scientific objectives.

| Window | What to do |
|---|---|
| *What data should be used for analysis?* | Abundance or biovolume, per taxonomic annotation, plankton group or trophic group. |
| *With what normalisation?* | None, log(x+1), double cube root, relative abundance, Hellinger, centred and reduced, or Box-Cox. |
| *Distance methods for PCoA?* | Euclidean (recommended after Hellinger transformation) or Bray-Curtis (recommended otherwise). |
| *Keep automatic dendrocut or select manually?* | Keep the automatic threshold of the hierarchical clustering, or enter your own. |

The module runs a PCA, a PCoA, a correspondence analysis and a hierarchical clustering, and projects the results on maps.

## 7. Where are the results?

All files are saved in the current MATLAB folder (the folder of your tsv files if you have just processed them):

| File or folder | Content |
|---|---|
| `base_<instrument>_<project>.mat` | Processed database (see [Outputs](outputs.md)) |
| `Map of stations.pdf` | Map of the samples |
| `raw analysis per station/` | One summary sheet per station |
| `global raw analysis/` | Bar charts and maps for the whole dataset (abundance or biovolume, diversity, dedicated taxa or groups) |
| `NBSS analysis/` | NBSS of all stations, truncated NBSS, map of NBSS slopes |
| `multivariates analysis/` | Outputs of the optional multivariate module |
| `<level>_excluded.xlsx` | Exclusion lists, to reload in a later run |

These figures are direct, unedited outputs, intended for data exploration. Check them and interpret them critically before using them in a publication.

## 8. Common problems

| Problem | Likely cause |
|---|---|
| *The folder contains no tsv files* | The export was not separated by sample, or the wrong folder was selected. |
| The pipeline stops while reading the files | A mandatory metadata field is missing or empty. Check the instrument page. |
| Abundances or biovolumes are implausible | A metadata field is in the wrong unit (e.g. mL instead of L). Check the instrument page. |
| A function is not found | The repository and `lib/` are not on the MATLAB path: run `addpath(genpath(pwd))` from the repository folder. |