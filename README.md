# aatomg - unlock the secrets of microproteins

This suite of programs performs reverse translation of peptide sequences to generate potential microgene nucleotide sequences. It then explores frameshift mutations in the coding genomic sequences, unveiling novel proteins. Ideal for researchers studying protein evolution, gene discovery, and functional genomics.

## Compilation

To compile the programs, navigate to the `microgene` directory and run `make`:

```bash
cd microgene
make
```

## Usage

```bash
./aatomg PATIFWY 1
```

- First argument: amino acid sequence as a string
- Second argument: strand flag
  - `0` = sense strand only
  - `1` = sense and antisense strands

**Note:** Make sure the data files `IsA.dat` and `Isma.dat` are present in the directory where you run `./aatomg`.

Zenodo DOI for the archived dataset: https://doi.org/10.5281/zenodo.18182315
