<p align="center">
  <img src="logo.png" alt="TS-SNAD logo" width="160">
</p>

<h1 align="center">TS-SNAD</h1>

<p align="center">
  <strong><code>TS-SNAD</code></strong><br>
  Tumor-Specific splicing neoAntigen detection
</p>

<p align="center">
  <em> A long-read transcriptomics workflow for discovering tumor-specific novel splice junctions and generating 9-mer neoantigens in nasopharyngeal carcinoma.</em>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.9%2B-3776ab" alt="Python 3.9+">
  <img src="https://img.shields.io/badge/pipeline-5_steps-0f766e" alt="5-step pipeline">
  <img src="https://img.shields.io/badge/focus-NPC_neojunctions-7c3aed" alt="NPC neojunctions">
</p>

<p align="center">
  <a href="#overview">Overview</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#script-reference">Script Reference</a> •
  <a href="#final-peptide-output-columns">Final Output</a> •
  <a href="#notes">Notes</a> •
  <a href="#04-pvactools">pVACbind</a>
</p>

---

`TS-SNAD` mainly consists of the following five standalone Python scripts:

1. quantify novel splice junction coverage across samples;
2. remove novel splice junctions that are present in GTEx normal tissues;
3. generate 9-mer peptides spanning novel splice junctions;
4. predict MHC-I binding affinity and filter candidate neoantigens;
5. check whether candidate peptides already exist in the reference human proteome.

The implementation tracks the workflow described in:

**Long-read sequencing reveals widespread novel splicing and neojunction-derived neoantigens in nasopharyngeal carcinoma**<br>
Yi Shuai, Hualiang Yao, Bo Wang, Grace TY Chung, Xiangeng Wang, Ming Zhong, Zhongxu Zhu, Cheuk Shuen Li, Chi Man Tsang, Kwok-Wai Lo*, Xin Wang*<br>
The Chinese University of Hong Kong

## Overview
Novel transcripts identified by long-read sequencing contain thousands of novel splice junctions (neojunctions), many of which are predominantly observed in tumor samples and are prioritized as tumor-specific. When these neojunctions occur within coding sequences (CDS), they can alter protein sequences and generate novel peptides that may be presented by HLA-I as potential neoantigens. TS-SNAD is designed to systematically identify and prioritize neojunction-derived neoantigens in five major steps.

```text
[01] identifying tumor-specific novel splice junctions using long-read transcriptomes with short-read support for splice junctions;
        |
        v
[02] removing novel splice junctions present in GTEx short-read and long-read splice junction matrices from normal tissues;
        |
        v
[03] performing CDS-aware translation and extracting splice junction-spanning 9-mer peptides;
        |
        v 
[04] predicting and filtering HLA-I binding affinity using pVACbind/NetMHCpan (e.g., IC50 < 500 nM and percentile rank < 0.5).
        |
        v
[05] checking whether predicted neoantigen peptides occur in the reference proteome and removing matched peptides.
        |
        v
neojunction-derived neoantigens
```

## Repository Layout

| Path | Role |
| --- | --- |
| `01_junction_coverage_tumor_specificity.py` | Builds a per-sample novel splice junction coverage matrix from STAR `SJ.out.tab` files and identifies tumor-specific novel splice junctions in the in-house cohort. |
| `02_gtex_junction_filter.py` | Screens novel splice junctions against GTEx short-read and long-read normal tissue splice junction references, excluding junctions detected in normal tissues.|
| `03_neoantigen_9mer_generator.py` | Generates novel splice junction-spanning 9-mer peptides from the CDS regions of novel isoforms.|
| `04_pvactools.py` | Merges the outputs of step 02 and step 03, then predicts MHC-I binding affinity using pVACbind/NetMHCpan and filters neojunction-derived neoantigens candidates based on binding thresholds.|
| `05_check_neoantigen_in_reference_proteome.py` | Optionally scans peptide sequences against a reference protein FASTA and labels peptides already present in the reference proteome. |
| `example/` | Bundled example inputs and example outputs for the three scripts. |
| `output/` | Default output folder used by the scripts. |

## Requirements

- Python 3.9 or newer
- `uv` recommended, or plain `pip`
- Core dependencies: `pandas`, `numpy`, `pyranges`, `biopython`, `pathos`

Install with `uv`:

```bash
cd TS-SNAD
uv sync
```

Or with `pip`:

```bash
pip install pandas numpy pyranges biopython pathos
```

## Script Reference

Before running the script commands below, enter the `TS-SNAD` directory:

```bash
cd TS-SNAD
```

### `01_junction_coverage_tumor_specificity.py`

This is the entry point for Step 1. It does two things in one pass:

1. Reads the long-read SQANTI3 junction file to extract novel splice junctions, and parses the STAR `SJ.out.tab` files for per-sample junction support. Then builds a per-sample junction coverage matrix;
2. quantifies tumor specificity and outputs the novel splice junctions that meet the criteria.

If you just want to run the bundled example, use defaults:

```bash
uv run python 01_junction_coverage_tumor_specificity.py
```

If you want to pass every path explicitly, there are two equivalent ways to provide sample groups:

1. Provide sample names directly on the command line:

```bash
uv run python 01_junction_coverage_tumor_specificity.py \
  --junction_file example/example_junctions.txt \
  --sj_path example/SJ.out.tab \
  --tumor_samples tumor1 tumor2 tumor3 tumor4 \
  --control_samples control1 control2 control3 \
  --output_dir ./output \
  --fold_change 10 \
  --max_normal_sum 5 \
  --min_tumor_sum 5
```

2. Provide `.txt` files containing sample names:

```bash
uv run python 01_junction_coverage_tumor_specificity.py \
  --junction_file example/example_junctions.txt \
  --sj_path example/SJ.out.tab \
  --tumor_samples example/tumor_samples.txt \
  --control_samples example/control_samples.txt \
  --output_dir ./output \
  --fold_change 10 \
  --max_normal_sum 5 \
  --min_tumor_sum 5
```

For a real cohort, the same command pattern looks like this after you replace the placeholders:

```text
uv run python 01_junction_coverage_tumor_specificity.py \
  --junction_file /path/to/junctions.txt \
  --sj_path /path/to/SJ.out.tab \
  --tumor_samples /path/to/tumor_samples.txt \
  --control_samples /path/to/control_samples.txt \
  --output_dir /path/to/results
```

Inputs:

| Argument | Default | What it expects |
| --- | --- | --- |
| `--junction_file` | `example/example_junctions.txt` | A SQANTI3 junction file containing splice junction category annotations. (`.txt`, tab-delimited). |
| `--sj_path` | `example/SJ.out.tab` | Directory containing one STAR `*.SJ.out.tab` file per sample. |
| `--tumor_samples` | `example/tumor_samples.txt` | Tumor sample names and/or existing `.txt` files containing sample names. | 
| `--control_samples` | `example/control_samples.txt` | Control sample names and/or existing `.txt` files containing sample names. |
| `--fold_change` | `10` | Multiplier used in `mean_tumor > fold_change * mean_normal`. |
| `--max_normal_sum` | `5` | Maximum allowed total novel splice junction coverage in controls. |
| `--min_tumor_sum` | `5` | Minimum required total novel splice junction coverage in tumors. |  

Notes for the input format:

- `--tumor_samples` and `--control_samples` support two input modes:
  - list sample names directly, e.g. `--tumor_samples tumor1 tumor2`
  - pass one or more `.txt` files, e.g. `--tumor_samples tumor_samples.txt`
- The same argument can mix names and `.txt` files, for example:
  - `--tumor_samples tumor1 tumor_more.txt`
  - `--control_samples control1 control_more.txt`
- `.txt` sample files are parsed by whitespace, so one-per-line or space-separated both work.
- Default sample groups are loaded from `example/tumor_samples.txt` and `example/control_samples.txt`.
- Coverage-matrix sample columns always follow this order: all tumor samples first, then all control samples.
- The script enforces a strict match between requested samples and `*.SJ.out.tab` files:
  - file count must equal requested sample count;
  - sample-name set must match exactly (no missing and no extra files).

The three threshold arguments above are applied together. A junction is kept only if all three conditions below are true:

```text
mean_tumor > fold_change * mean_normal
AND total_normal_reads < max_normal_sum
AND total_tumor_reads  > min_tumor_sum
```

Outputs written to `--output_dir`:

| File | Meaning |
| --- | --- |
| `novel_junctions_matrix.csv` | Novel splice junction coverage matrix with one row per junction and one column per sample. This file is standard comma-delimited CSV. |
| `tumor_specific_novel_junctions.csv` | Tumor-specific novel splice junctions. This file is tab-delimited even though the extension is `.csv`. |

With the example command above, the script writes:

- `output/novel_junctions_matrix.csv`
- `output/tumor_specific_novel_junctions.csv`

The output directory is created automatically if it does not already exist.

### `02_gtex_junction_filter.py`

This script checks the tumor-specific junction list from Step 1 against two GTEx references:

- GTEx short-read junctions from the GCT `Name` column
- GTEx long-read junctions from the FLAIR GTF coordinates

It writes two "junction exists in GTEx" reports, then writes a third file containing the tumor-specific junctions that were not found in either GTEx dataset.

Default example run:

```bash
cd TS-SNAD
tar -xzf example/GTEx/flair_filter_transcripts.gtf.tgz
uv run python 02_gtex_junction_filter.py
```

The GTEx short-read GCT is a large file; a download link will be added in a future update. Place it at `example/GTEx/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct` before running Step 2. The long-read GTEx example annotation is tracked as `example/GTEx/flair_filter_transcripts.gtf.tgz` to keep the repository smaller. Extract it before running Step 2; pass the decompressed `.gtf` file to `--gtex_lr_file`, not the `.tgz` archive.

Custom run:

```bash
uv run python 02_gtex_junction_filter.py \
  --gtex_sr_file example/GTEx/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct \
  --gtex_lr_file example/GTEx/flair_filter_transcripts.gtf \
  --junction_file example/example_junctions.txt \
  --tumor_specific_junctions output/tumor_specific_novel_junctions.csv \
  --output_dir ./output
```

Inputs:

| Argument | Default | What it expects |
| --- | --- | --- |
| `--gtex_sr_file` | `example/GTEx/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct` | GTEx short-read splice-junction matrix in GCT format. |
| `--gtex_lr_file` | `example/GTEx/flair_filter_transcripts.gtf` | GTEx long-read FLAIR transcript annotation in GTF format. |
| `--junction_file` | `example/example_junctions.txt` | SQANTI3 junction file (same as Step 1). |
| `--tumor_specific_junctions` | `output/tumor_specific_novel_junctions.csv` | Tumor-specific novel splice junctions (output from Step 1).  |
| `--output_dir` | `./output` | Folder for the GTEx reports and filtered list. |

Outputs written to `--output_dir`:

| File | Meaning |
| --- | --- |
| `GTEx_cohort_exist.csv` | Tumor-specific novel splice junctions that can be found in the GTEx short-read splice junction matrix. |
| `GTEx_LR_cohort_exist.csv` | Tumor-specific novel splice junctions  that can be found in the  GTEx long-read splice junction matrix. |
| `{tumor_specific_junctions_basename}_gtex_filter.csv` | Tumor-specific novel splice junctions after removing any junction found in either GTEx dataset. | 

Examples of the third output name:

- `tumor_specific_novel_junctions_gtex_filter.csv`

With the default arguments, the script writes:

- `output/GTEx_cohort_exist.csv`
- `output/GTEx_LR_cohort_exist.csv`
- `output/tumor_specific_novel_junctions_gtex_filter.csv`

### `03_neoantigen_9mer_generator.py`

This script takes transcript annotation plus sequence data and generates every 9-mer peptide that spans a novel splice junction.

For each transcript, it:

1. reads exon and CDS intervals,
2. removes intronic positions,
3. translates the coding sequence,
4. maps junction boundaries into amino-acid space,
5. extracts the sliding 9-mer windows that cross the junction.

Default example run:

```bash
uv run python 03_neoantigen_9mer_generator.py
```

Custom run:

```bash
uv run python 03_neoantigen_9mer_generator.py \
  --cds_gff example/example.cds.gtf \
  --junction_file example/example_junctions.txt \
  --dna_seq example/example.fasta \
  --protein_seq example/example.faa \
  --output_file ./output/positive_novel_junction_peptides_example.csv
```

Parallel mode:

```bash
uv run python 03_neoantigen_9mer_generator.py \
  --cds_gff example/example.cds.gtf \
  --junction_file example/example_junctions.txt \
  --dna_seq example/example.fasta \
  --protein_seq example/example.faa \
  --output_file ./output/positive_novel_junction_peptides_example.csv \
  --parallel
```

Inputs:

| Argument | Default | What it expects |
| --- | --- | --- |
| `--cds_gff` | `example/example.cds.gtf` | Transcript exon/CDS annotation. The bundled example is GTF-formatted. |
| `--junction_file` | `example/example_junctions.txt` | SQANTI3 junction file (same as Step 1). |
| `--dna_seq` | `example/example.fasta` | Transcript DNA FASTA. |
| `--protein_seq` | `example/example.faa` | Transcript protein FASTA. |
| `--output_file` | `./output/positive_novel_junction_peptides_example.csv` | Output 9-mer peptide table. |
| `--parallel` | off | Enables parallel processing with `pathos`. |

Outputs:

| File | Meaning |
| --- | --- |
| `positive_novel_junction_peptides_example.csv` | Tab-delimited table of 9-mer peptides spanning novel junctions. |

With the default arguments, the script writes:

- `output/positive_novel_junction_peptides_example.csv`
The Step 3 output file is tab-delimited. These columns are the ones you will usually inspect first:

| Column | Meaning |
| --- | --- |
| `transcript_id` | Isoform identifier from the input FASTA/GTF/GFF. |
| `new_peptide` | A 9-mer spanning the junction, or a sentinel value such as `junc-loc-in-start-or-end-code`. |
| `junction_id` | Numeric junction identifier within the isoform. |
| `strand` | Transcript strand. |
| `start_acid`, `end_acid` | Amino-acid window boundaries used to derive the peptide set. |
| `start_to_end` | Amino-acid segment around the junction before 9-mer sliding is applied. |
| `acid` | Full translated coding peptide used for validation. |
| `mrna_seq` | CDS nucleotide sequence used for translation. |
| `isoform_seq` | Full isoform nucleotide sequence. |
| `exon_flag` | Encoded exon segments in `id:start-end;` form. |
| `cds_flag` | Encoded CDS-covered segments in `id:start-end;` form. |
| `junc_flag_for_acid` | Junction-aware flag positions in amino-acid space. |



### `04_pvactools.py`

This script takes the outputs of Step 2 and Step 3, intersects them to retain only tumor-specific novel splice junction-derived 9-mers, runs pVACbind to predict MHC-I binding affinity.

It does four things in sequence:

1. intersects the GTEx-filtered tumor-specifc novel junction list (output from Step 2) with the 9-mer peptide table (output from Step 3) by `(isoform, junction_id)`; saves the surviving rows as `neoantigen.intersected.tsv` and writes the surviving unique 9-mers as `neoantigen.pep.fa`;
2. runs `pvacbind run` with NetMHCpan to predict MHC-I binding affinity;
3. runs `pvacbind binding_filter` to keep only strong binders (IC50 < 500 nM, median percentile < 0.5);
4. matches `MHC_Class_I/<prefix>.MHC_I.all_epitopes.filtered.tsv` column `Epitope Seq` back to `neoantigen.intersected.tsv` column `new_peptide`, then saves the matched full rows as `neoantigen.intersected.filtered_epitopes.tsv`.

Custom run:

```text
rm -rf ./pvactools_output
uv run python 04_pvactools.py \
  --tumor_specific_junctions output/tumor_specific_novel_junctions_gtex_filter.csv \
  --neoantigen_file output/positive_novel_junction_peptides_example.csv \
  --hla_alleles example/hla_alleles.txt \
  --output_directory ./pvactools_output \
  --prefix NPC
```

`--hla_alleles` supports either:

- a comma-separated list such as `HLA-B*44:02,HLA-B*44:03`
- a text file path with one HLA allele per line

Inputs:

| Argument | Required | What it expects |
| --- | --- | --- |
| `--tumor_specific_junctions` | yes | Step-02 GTEx-filtered junction CSV (columns: `isoform`, `junction_num`). |
| `--neoantigen_file` | yes | Step-03 9-mer neoantigen TSV (columns: `transcript_id`, `junction_id`, `new_peptide`, ...). |
| `--hla_alleles` | yes | Comma-separated HLA-I alleles or a txt file path with one HLA allele per line, e.g. `HLA-A*02:07,HLA-B*46:01,HLA-C*01:02` or `example/hla_alleles.txt`. |
| `--output_directory` | yes | Root output directory; pVACbind writes into `<output_directory>/MHC_Class_I/`. |
| `--prefix` | no | Sample prefix for pVACbind output naming (default: `NPC`). |

Outputs written to `--output_directory`:

| File | Meaning |
| --- | --- |
| `neoantigen.intersected.tsv` | The full Step-03 rows retained after the step-02 × step-03 intersection, before 9-mer deduplication. |
| `neoantigen.intersected.filtered_epitopes.tsv` | Rows from `neoantigen.intersected.tsv` whose `new_peptide` exactly matches a filtered pVACbind `Epitope Seq`. |
| `neoantigen.pep.fa` | FASTA of unique valid 9-mer peptides that survived the step-02 × step-03 intersection. |
| `pvactools.log` | Full stdout/stderr log from `pvacbind run`. |
| `MHC_Class_I/<prefix>.MHC_I.all_epitopes.tsv` | pVACbind raw binding affinity predictions for all alleles. |
| `MHC_Class_I/<prefix>.MHC_I.all_epitopes.filtered.tsv` | Predictions after binding filter (IC50 < 500 nM, percentile < 0.5). |
| `MHC_Class_I/binding_filter.log` | Log from `pvacbind binding_filter`. |

Notes for Step 4:

- Requires `pvacbind` to be installed in the current environment.
- HLA alleles must use the standard colon format: `HLA-A*02:07`, not `HLA-A0207`.
- Only 9-mers composed of the 20 standard amino acids are kept; error-tag rows (e.g. `junc-loc-in-start-or-end-code`, `slide-error`) are discarded automatically.

### `05_check_neoantigen_in_reference_proteome.py`

Run Step 5 after Step 4:

```bash
uv run python 05_check_neoantigen_in_reference_proteome.py \
  pvactools_output/neoantigen.intersected.filtered_epitopes.tsv \
  example/uniprotkb_Human_AND_model_organism_9606_2024_07_25_sp.fasta \
  --peptide-column-name new_peptide \
  --output-prefix output/neoantigen_reference_check
```

The example command uses the UniProt human reference proteome FASTA at `example/uniprotkb_Human_AND_model_organism_9606_2024_07_25_sp.fasta`. Replace this path if you want to screen against a different reference proteome.

## Notes

- Step 1 and Step 2 are chained by file: Step 2 reads the tumor-specific output from Step 1.
- Step 3 does not read the Step 2 filtered file directly. It generates 9-mers for all novel junctions in the SQANTI3 table. Step 4 performs the intersection.
- Step 4 chains Step 2 and Step 3: it takes both outputs and retains only 9-mers whose junction survived the GTEx filter before running pVACbind.
- `02_gtex_junction_filter.py` and `03_neoantigen_9mer_generator.py` write tab-delimited content with a `.csv` extension. Keep that in mind if you load them outside pandas.
- `03_neoantigen_9mer_generator.py --parallel` creates a `pathos` process pool with a fixed size of 150 workers. That is fine on a large server and a bad default on a laptop.
- `04_pvactools.py` now executes `pvacbind` directly from the current environment. If you use `uv`, prefer `uv run python 04_pvactools.py ...`.
