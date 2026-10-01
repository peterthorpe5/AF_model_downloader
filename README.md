# Standalone AlphaFold Database downloader

This Python 3.10+ command retrieves **published AlphaFold Database models**;
it does not predict structures. It maps protein accessions through UniProt
before requesting a canonical AlphaFold Database F1 model. It needs only the
Python standard library.

Supported `--id-type` choices:

| Type | Example | Input |
| --- | --- | --- |
| `ncbi-protein` | `AAB61673.1`, `NP_009225.1` | A mixed list of INSDC GenBank/ENA/DDBJ and RefSeq **protein** accessions |
| `genbank-protein` | `AAB61673.1` | INSDC GenBank/ENA/DDBJ CDS **protein** accessions |
| `refseq-protein` | `NP_000537.3` | RefSeq **protein** accessions |
| `gene-id` | `7157` | Numeric NCBI Gene IDs; one gene may yield several UniProt proteins |
| `uniprot` | `P04637` | UniProt accessions without an ID mapping job |

For `ncbi-protein`, the command automatically submits separate UniProt mapping
jobs for RefSeq proteins and INSDC protein accessions. It retains **every**
UniProt match and records both accession types in `manifest.tsv`. A valid
mapping does not prove that the original protein has the same sequence as the
UniProt model. Use a protein FASTA to check sequence identity before download.

## Quick start: one or more accessions

Unpack the ZIP, or git clone https://github.com/peterthorpe5/AF_model_downloader.git  and change to the `AF_model_downloader` directory.

```bash
python3 fetch_alphafold_models.py \
  --id-type ncbi-protein \
  --id AAB61673.1 \
  --id NP_009225.1 \
  --output-dir "$HOME/alphafold_mixed_example" \
  --include-pae
```

For a list, make a UTF-8 text file with **one accession per line and no
header**. Blank lines and lines beginning with `#` are ignored:

```text
AAB61673.1
NP_009225.1
```

```bash
python3 fetch_alphafold_models.py \
  --id-type ncbi-protein \
  --ids-file "$HOME/my_protein_accessions.txt" \
  --output-dir "$HOME/alphafold_mixed_batch" \
  --include-pae
```

Output directories must be **new or empty**. The command refuses to overwrite
an earlier run. The default limit is 50 distinct models; use `--max-models 100`
to raise it when needed (maximum 500). For a different request, pick a new
output directory.

## Preferred: supply the actual protein sequences

If you have a multi-FASTA, you can use it as the **only input file**. Each
header must start with the NCBI **protein** accession, preferably with its version;
the remaining description is ignored. Sequences must be unaligned. Use the
real sequences from NCBI or your analysed data, not the abbreviated example
below:

```text
>AAB61673.1 protein description
THE_COMPLETE_PROTEIN_SEQUENCE_HERE
>NP_009225.1 protein description
THE_COMPLETE_PROTEIN_SEQUENCE_HERE
```

```bash
python3 fetch_alphafold_models.py \
  --id-type ncbi-protein \
  --protein-fasta "$HOME/my_proteins.faa" \
  --output-dir "$HOME/alphafold_checked_batch" \
  --include-pae
```

The command compares each submitted sequence with the full UniProt sequence
reported by AlphaFold Database. It downloads a model for an input only when
the sequences match **exactly**. If they differ, the manifest keeps the ID
conversion but reports `SEQUENCE_MISMATCH` and writes no input-named model for
that record. This is a mapping check, not a coordinate-level residue check.
For an accession list and a separate FASTA, use `--ids-file` together with
`--sequence-fasta`; the FASTA must contain exactly the same accession headers.

**Ken's example needs care:** UniProt identifies `AAB61673.1` as a BRCA1
cross-reference but annotates its translated sequence as different from the
curated BRCA1 sequence, with an erroneous CDS choice. A mapping to `P38398`
must therefore **not** be described as an exact structure for `AAB61673.1`
without checking the actual sequences. The FASTA mode will flag a mismatch
against the canonical model rather than silently substituting it. With an
accession-only list, the manifest says `NOT_PROVIDED` for the sequence check.

There is no guarantee that an accession has a UniProt mapping or that a
mapped UniProt protein has a published canonical model. The script reports
these separately as `UNMAPPED` and `MODEL_NOT_AVAILABLE`; it does not treat a
related protein as a verified structural match.

## Files and naming

```text
OUTPUT_DIR/
  manifest.tsv
  models/
    AF-<UniProt>-F1-model_v<version>.cif
    AF-<UniProt>-F1-predicted_aligned_error_v<version>.json
  by_input/
    <submitted ID>__<UniProt>__AF-<UniProt>-F1-model_v<version>.cif
    <submitted ID>__<UniProt>__AF-<UniProt>-F1-predicted_aligned_error_v<version>.json
```

`manifest.tsv` is the accession conversion table. It records `input_id`,
`mapping_source`, `uniprot_accession`, `status`, the canonical `model_path`,
the accession-labelled `named_model_path`, the corresponding PAE paths when
requested, SHA-256 checksums, model version, sequence check and error details.
Each mapped UniProt accession gets its own row. The `by_input/` names preserve
the accession supplied by Ken or Maddy and identify the actual UniProt model.
The canonical file is downloaded once per UniProt accession; input-named
files use hard links where possible or a copy when links are unavailable.
The model path stays empty for unmapped, unavailable or mismatching input IDs.

Other statuses: `DOWNLOADED` means a model was written; `FAILED` means an API,
validation or file error. Inspect `detail` and rerun into a fresh directory.
The command exits `0` when all requests have a known outcome (including
unmapped or unavailable models), `1` for failed requests or requested PAE
failures, and `2` for invalid input or another fatal error. Use `--help` for
all available options.

These are AlphaFold Database predictions with pLDDT and, where available,
PAE. The downloads do not have AlphaFold 3 pTM, ipTM or ranking scores.

## Test

```bash
python3 -m unittest discover -s tests -v
```

The included tests use controlled HTTP responses and make **no network
calls**. The authoring environment could not make a live UniProt/AlphaFold
request. Run a small smoke test on a machine with internet access and inspect
`manifest.tsv` before relying on any mapping.

Sources: [UniProt's current ID mapping fields](https://rest.uniprot.org/configure/idmapping/fields),
[NCBI accession formats](https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/),
[UniProt's annotation of AAB61673.1](https://www.uniprot.org/uniprotkb/P38398/entry),
and the [AlphaFold Database API description](https://academic.oup.com/nar/article/50/D1/D439/6430488).
