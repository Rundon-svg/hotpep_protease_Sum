# Peptidase Annotation Summary Tool

This Python script summarizes peptidase annotation results from multiple species using data from the **Hotpep_protease** tool for my project.

## Input Data Structure

The test input data comes from Hotpep_protease(https://www.sciencedirect.com/science/article/pii/S2666952820300431), with each species' protein sequence files (longest deduplicated collections via CD-HIT) resulting in one output folder per species, obtained using Hotpep's default parameters (except viral).

The input consists of folders for this script for each species, each containing a `peptidases` directory with a `summary.txt` file. The `summary.txt` file must include the following columns:

- **Merops family**: Classification of the peptidase.
- **proteins**: Count of proteins associated with each Merops family.

The folder structure should be as follows:

```
input_path/
    species1/
        peptidases/
            summary.txt
    species2/
        peptidases/
            summary.txt
    ...
```

## Usage

To run the script, use the following command:

```bash
python summarize_peptidase_results.py -in <input_path>
```

Replace `<input_path>` with the path to the folder containing all species results.

## Output

The script generates two output files:

1. **combined_summary.csv**: Gene counts of each peptidase families/subfamilies for all species.
2. **summary_statistics.csv**: Total counts categorized by Merops family and subfamily, ensuring all categories are represented.
   Statistics by three categories:  
   - **Total**
   - **9 catalytic types**: Aspartic (A) Peptidases, Cysteine (C) Peptidases, Glutamic (G) Peptidases, Metallo (M) Peptidases, Asparagine (N) Peptide Lyases, Mixed (P) Peptidases, Serine (S) Peptidases, Threonine (T) Peptidases, Peptidases of Unknown (U) Catalytic Type
   - **Family**

## Requirements

- Python 3.x
- Pandas library

Install the required library using:

```bash
pip install pandas
```

## Usage


```bash
python hotpep_statistics_multi_en.py -in /path/to/input_folder
```

This will process the data and output the results in the same directory.
