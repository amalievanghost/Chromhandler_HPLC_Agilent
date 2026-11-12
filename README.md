# Chromhandler_HPLC_Agilent

Preprocess Agilent HPLC peak reports (`.txt`) into concentration time courses, saved as [EnzymeML](https://enzymeml.org) documents and plotted per replicate.

Built on [chromhandler](https://github.com/haeussma/chromhandler) and [pyenzyme](https://github.com/EnzymeML/PyEnzyme). Currently set up for cGAS reactions (ATP + GTP → cGAMP), but any molecule with a calibration can be added.

## Workflow

Quick start

For a new experiment:

1) Add reaction conditions in the metadata table under: Input/metadata.xlsx
2) Sort the raw files run: Sorting.ipynb
3) Calibrate with Calibration.ipynb (only if the standards or molecules changed)
4) Analyse with Analysis.ipynb → concentrations + plots will be saved in Results/

Look for `# CHANGE HERE #` in the notebooks: those are the lines you edit for a new experiment or a new molecule.
There are currently two different styles for the plot with broken axis, plot styles can be changed in the .py scripts under Input/

### 1. Add the reaction conditions to `Input/metadata.xlsx`

One row per **condition**, not per replicate. Replicates `_1`, `_2`, `_3` automatically get the metadata of their condition.

| experiment_name | reaction_vol [mL] | temperature [C] | pH | HEPES [mM] | TrisHCl [mM] | MgCl2 [mM] | tmmcGAS [mg/ml] | thscGAS [mg/ml] | ATP [mM] | GTP [mM] | … |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20260520_100mM-HEPES_0.04_tmmcGAS | 1 | 37 | 7.0 | 100 | | 10 | 0.04 | | 0.5 | 0.5 | |

>[!IMPORTANT]
>Formatting and guidelines:
>- `experiment_name` is the raw-data folder name, **without** `_1/_2/_3`.
>- Units go in square brackets in the header. They are stripped when the table is read.
>- Column names (without unit) must match the molecule/protein ids **exactly**, including upper/lower case (`TrisHCl`, not >`TrisHCL`).
>- Empty cell = 0. Leave molecules you didn't use empty.
>- Small molecules are in **mM**; proteins in **mg/mL**. Protein concentrations are converted to µM in `Analysis.ipynb` >using the molecular weights in `MW`. Add the MW there if you add a new protein.

### 2. Sort the raw data: `Sorting.ipynb`

Put the exported Agilent reports of one experiment in `Raw-data/<experiment_name>/`.

>[!IMPORTANT] **Naming**
>
>```
>Raw-data/
>└── YYYYMMDD_conditions/                      # folder = experiment_name
>    ├── YYYYMMDD_1_conditions_0_min.txt       # replicate 1
>    ├── YYYYMMDD_1_conditions_10_min.txt
>    ├── ...
>    ├── YYYYMMDD_2_conditions_0_min.txt       # replicate 2
>    ├── YYYYMMDD_3_conditions_0_min.txt       # replicate 3
>    ├── ...
>    └── YYYYMMDD_NEG_conditions_1440_min.txt  # negative control
>```
>
>Example: `20260520_1_100mM-HEPES_0.04_tmmcGAS_10_min.txt`
>
>- `YYYYMMDD`: date of the experiment. Optional in the file names, but always 8 digits if present. In the folder name, >YYYYMMDD makes experiments sort chronologically.
>- `1`, `2`, `3`: replicate number. It decides which replicate folder the file is sorted into.
>- `NEG`: negative control. It is not sorted and stays in `Raw-data`.
>- `conditions`: free text describing the reaction, e.g. `100mM-HEPES_0.04_tmmcGAS`. It is only a label; the values used in >the analysis come from `metadata.xlsx`.
>- `<time>_min`: reaction time at sampling. chromhandler reads the time from here, so it must be the **last number before >`_min`**.

Set `SOURCE_PATH` to that folder and run the notebook. It copies the files into

```
Data/time_sample/<experiment_name>/
├── <experiment_name>_1/
├── <experiment_name>_2/
└── <experiment_name>_3/
```

Files without a replicate number, such as the negative controls (`NEG`), stay in `Raw-data`.

### 3. Calibrate: `Calibration.ipynb`

Defines every small molecule and protein once and saves it as JSON:

- **Proteins** (UniProt ID) → `Results/Proteins/`
- **Constant molecules** such as buffer and salts (PubChem CID, no calibration) → `Results/Molecules/`
- **Calibrated molecules** (substrates and products) → linear standard curve from the standards in `20260610_Standards/<molecule>/` → `Results/Molecules/`

To add a molecule:
- **Constant:** define it in the constant-molecules cell and add it to the save list.
- **Calibrated:** add a standards folder named after the molecule, and add it to `target_RT`, `tolerances` and `PUBCHEM_CIDS`.

In both cases, also add a matching column to `metadata.xlsx`.

Calibration only needs to be rerun when the standards or the species change.

### 4. Analyse: `Analysis.ipynb`

1. Loads the molecules and proteins from `Results/`. Retention times can be fine-tuned per molecule here.
2. Loads `Input/metadata.xlsx` and expands it to one row per replicate.
3. Set `sample_dir` to `Data/time_sample/<experiment_name>`. For each replicate, it reads the chromatograms, adds the initial concentrations and assigns the peaks.
4. Converts peak areas to concentrations and saves an EnzymeML document to `Results/enzymeml/<experiment_name>.json`.
5. Plots all replicates and mean ± SD with a broken y-axis and saves the PNGs to `Results/`. Change the plot layout in `Input/broken_open_frame_axis.py`.

If a molecule doesn't change during the reaction (buffer, salts), add it to the constant list in the analysis loop.

## Folder structure

```
├── Input/
│   ├── metadata.xlsx                 # reaction conditions, one row per condition
│   └── broken_*_frame_axis.py        # plot layout
├── Raw-data/<experiment_name>/       # Agilent .txt reports as exported
├── Data/time_sample/<experiment_name>/<experiment_name>_1..3/   # sorted replicates
├── 20260610_Standards/<molecule>/    # calibration standards
├── Results/
│   ├── Molecules/  Proteins/         # species JSON files from calibration
│   ├── enzymeml/                     # EnzymeML documents
│   └── *.png                         # plots
├── utility/                          # plotting helpers
├── Sorting.ipynb
├── Calibration-JSON.ipynb
└── Analysis.ipynb
```

## Installation

```bash
conda create -n chromhandler python=3.11
conda activate chromhandler
pip install chromhandler pyenzyme pandas numpy matplotlib brokenaxes openpyxl jupyter
```
