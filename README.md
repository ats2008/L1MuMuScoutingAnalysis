[README UPATED WITH Google-Antigravity, please pardon the errs. !] 
# L1 Scouting Dimuon Analysis

This repository contains the analysis framework for L1 Scouting Dimuon physics. It performs event selection, categorization, and histogramming using Python-based tools centered around `uproot`, `awkward`, and `vector`.

## Table of Contents
- [Project Structure](#project-structure)
- [Prerequisites](#prerequisites)
- [Workflow Overview](#workflow-overview)
- [Running the Analysis](#running-the-analysis)
    - [Configuration](#configuration)
    - [Step 1: Histogram Production](#step-1-histogram-production)
    - [Step 2: Merging Outputs](#step-2-merging-outputs)
    - [Step 3: Summary and Plotting](#step-3-summary-and-plotting)
- [Detailed Script Reference](#detailed-script-reference)

---

## Project Structure

```text
.
├── data/                   # Configuration files and input file lists
│   ├── v1_files.json       # Dictionary mapping process names to input ROOT files
│   └── xsdb.json           # Cross-section database
├── python/                 # Source code for analysis and utilities
│   ├── analysis.py         # Core event loop and histogramming logic
│   ├── runHistMaker.py     # Driver script to submit/run histogram jobs
│   ├── runAnalysis.py      # Driver script to analyze summaries and yields
│   ├── util.py             # Library for physics selections (Isolation, DeltaR, etc.)
│   ├── mergeHistos.py      # Script to merge JSON outputs
│   ├── printSummary.py     # Script to print yields and efficiencies
│   └── ...                 # Other plotting and export scripts
└── README.md               # This documentation
```

## Prerequisites

The analysis runs on **Python 3**. The following packages are required:

- `uproot`: For reading ROOT files.
- `awkward`: For manipulation of jagged arrays and columnar data.
- `vector`: For 4-vector kinematics.
- `numpy`: For numerical operations.
- `tabulate`: For printing summary tables.

You can install them via pip:
```bash
pip install uproot awkward vector numpy tabulate
```

## Workflow Overview

The standard analysis chain is as follows:
1.  **Define Inputs**: List your source ROOT files in a JSON configuration (e.g., `data/v1_files.json`).
2.  **Run Analysis**: Use `runHistMaker.py` to run `analysis.py` over the input files. This produces JSON files containing histograms (one per input file/chunk).
3.  **Merge Results**: Use `mergeHistos.py` to combine the scattered JSON files into a single output per process.
4.  **Analyze Results**: Use `runAnalysis.py` (or `printSummary.py`) to check yields, efficiencies, and produce plots.

---

## Running the Analysis

### Configuration
Update `data/v1_files.json` with your dataset paths. The format should be:
```json
{
    "ProcessName_1": ["/path/to/file1.root", "/path/to/file2.root"],
    "ProcessName_2": ["/path/to/fileA.root"]
}
```

### Step 1: Histogram Production
The `runHistMaker.py` script is the main driver for producing histograms. It reads the file list and generates commands to run `analysis.py` for each file.

**Dry Run (Print commands only):**
```bash
python3 python/runHistMaker.py --ver v1
```

**Execute Analysis:**
Use the `--exec` flag to actually run the commands.
```bash
python3 python/runHistMaker.py --ver v1 --exec
```

**Common Options:**
- `--ver <version>`: Specify a version tag (creates output in `results/analysis/<version>/`).
- `--procs <ProcessName>`: Run only for a specific process (e.g., `MinBias`).
- `--isTest`: Run only on a small subset of events (for debugging).

### Step 2: Merging Outputs
After running the analysis, outputs are typically split across multiple JSON files. Merge them using `mergeHistos.py`.

```bash
python3 python/mergeHistos.py \
    -s "results/analysis/v1/MinBias/parts/out_data_highptMuMu_MinBias_*.json" \
    -o "results/analysis/v1/MinBias/out_data_highptMuMu_MinBias.json"
```

*Note: The script supports replacing `@` with `*` in the search string for convenience.*

### Step 3: Summary and Plotting
To view the cutflow summary (Yields and Efficiencies), use `runAnalysis.py` (which runs the summary analyzer) or `printSummary.py`.

```bash
python3 python/printSummary.py \
    -i results/analysis/v1/MinBias/out_data_highptMuMu_MinBias.json \
    -p MinBias
```

---

## Detailed Script Reference

### `python/analysis.py`
The core worker script. It processes a single ROOT file and produces a JSON file with histograms.
- **Inputs**:
    - `-i, --inputFile`: Path to the input ROOT file.
    - `-t, --tag`: Tag/Label for the dataset.
    - `-d, --destination`: Output folder.
    - `-e, --eMax`: Max events to process (default: all).
- **Core Logic**:
    - Selects events based on Dimuon criteria (defined in `util.py`).
    - Calculates isolation variables (Puppi, PF).
    - Can categorize events into `inclusive`, `highptMuMu`, etc.
    - Saves histograms dictionary to JSON.

### `python/runHistMaker.py`
Driver script to manage batched histogram production.
- **Functionality**:
    - Reads input file lists from `data/v1_files.json`.
    - Iterates over specified processes and generates specific `analysis.py` commands for each file chunk.
    - Constructs output paths based on provisioned version (`--ver`).
- **Inputs**:
    - `-v, --ver`: Analysis version (controls output path).
    - `--exec`: If not present, only prints the commands.
    - `-p, --procs`: Filter for specific process names in the config file.
    - `--isTest`: Run on a small subset for debugging.

### `python/runAnalysis.py`
Script to analyze aggregated results and produce summaries.
- **Functionality**:
    - Aggregates the scattered JSON results produced by `runHistMaker.py`.
    - Iterates over defined categories (e.g., `highptMuMu`, `inclusive`).
    - Calls `printSummary.py` internally to extract yields and efficiencies.
    - Generates a consolidated `analysis_summary.json` file.
    - Prints formatted tables (Success/Yield/Efficiency) to the console.
- **Inputs**:
    - `--exec`: Execute the internal summary commands.
    - `-o, --ofname`: Output summary filename.
    - `--tag`: Tag for the output summary.

### `python/util.py`
Contains the physics logic and selection definitions.
- **Key Functions**:
    - `addSelectedDimuons`: Applies kinematic cuts, ID, and isolation.
    - `addGenSelectedDimuons`: Selection for Gen-matched objects.
    - `addPuppiIsolation` / `addPFIsolations`: Computes isolation sums in cones.
    - `getHistograms`: Helper to format data for JSON export.

### `python/mergeHistos.py`
Merges JSON histogram files.
- **Inputs**:
    - `-s, --searchString`: Glob pattern for input files.
    - `-o, --outFile`: Path for the merged output.

### `python/printSummary.py`
Prints a tabulated summary of the analysis results.
- **Inputs**:
    - `-i, --inputFile`: Merged JSON file.
    - `-p, --proc`: Process name (for labeling).
    - `--printAllHistograms`: Dumps details of all histograms found.
