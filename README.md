# oligocompose

##### Jamie Heather | CCR @ MGH | 2024

[![License](https://img.shields.io/github/license/JamieHeather/oligocompose?label=license)](./LICENSE)

The aim of this repo is to facilitate the production of large oligo pool library sequences, particularly potential MHC peptide epitope libraries, from shorter substrings.

* Information of sequences to be generated should be supplied in a tab-separated file.
  * The first column must contain a name or identifier (which needn't be unique).
  * All additional columns contain the sequences to combine in the output oligos, with optional additional comma-delimited arguments providing additional information about each sequence in square brackets after the sequence.
  * Different lines may have different numbers of columns; the input file needn't be a 'proper' whole-grid tab-separated file. 
  * All lines must be left aligned, i.e. with no empty columns separating sequences. 
  * Comment lines that will not contribute to the output can also be included by beginning the line with a `#` character.
  * This file is referred to via the `-in` or `--in_file` flags, e.g.:

```bash
oligocompose -in my_input_file.tsv
```


* Sequences can be provided via one of several other ways, in any combination within a line:
  * Explicitly as a simple DNA sequence, or as a amino acid sequence, which will be reverse-translated using the most commonly used codon to encode that residue in humans.
    * While the script attempts to automatically infer the sequence type (defaulting to nucleotides), it's preferable to explicitly label each sequence type in the square bracket parameter fields with `[nt]` or `[aa]`.
    * E.g. `TTT	aaa[nt]	CCC` produces 'TTTAAACCC', while `TTT	aaa[aa]	CCC` produces 'TTTGCCGCCGCCAAA', with the codons for three alanines encoded between the TTT and CCC sequences.


* Individual sequences that are to be used frequently (e.g. conserved adapter sequences) can be included in a reference 'fixed sequence' file, which can be specified using the `-f` or `--fixed_file` flag. 
  * This should fit the criteria of the input oligo tsv file, except it should be only two columns, name and reference sequence.
  * Those sequences can then be referred to in the input oligo file but using those names, bracketed with ampersands (e.g. `&name_of_sequence&`).

```bash
oligocompose -in my_input_file.tsv -f fixed_sequences.tsv
```

* For variable sequences (e.g. different sequences from a pool), instead of explicitly providing each individual a DNA string, users can provide a path to a file containing a list of sequences can be used, bracketed with dollar signs (e.g. `$path_to/some_file.tsv$`)
  * This file must be formatted the same as the fixed sequence file, i.e. a two column tsv of `name\tsequence`.
  * In this situation an output sequence will be produced for every sequence in the additional file. E.g. when `some_file` contains the sequences AAA, CCC, GGG, the line `TT    $some_file$    TTT` would produce TTAAATTT, TTCCCTTT, and TTGGGTTT.
  * Note that this file cannot contain links to other files, only explicit sequences or references to them.


* Degenerate DNA bases (using IUPAC codes) can also be used.
  * E.g. `TT NNN AA` would produce 64 oligos, with each combination of all four nucleotides at each of the three degenerate positions.


* Instead of just inserting whole sequences (provided either directly or via a file), all sub-sequences of a specified length can be individually combined. 
  * This is particularly useful for generating oligos tiled along the length of a gene or protein sequence.
  * After giving a sequence, identifier, or path, integer values of length of subsequence to take can be specified in square brackets.
  * Note that this behaviour defaults to a sliding window of 1 for amino acids, and 3 for nucleotides, assuming that its been provided with an in-frame ORF and that codons are to be maintained. This behaviour can be altered by providing any integer to the script via the `--sl / --step_length` flag.
  * E.g. `ACDEFGHIKLMNPQRSTVWY[10,aa]` would produce nucleotides encoding every possible 10-mer across this theoretical protein sequence (`ACDEFGHIKL`, ` CDEFGHIKLM`, `DEFGHIKLMN`, ... `MNPQRSTVWY`)
    * Running the same code with `-sl 5` included in the script would instead produce 10-mers every 5 residues apart (just `ACDEFGHIKL`, ` GHIKLMNPQR`, and `MNPQRSTVWY`).
  * Alternatively if running on nucleotide sequences, coding sequences can be maintained by providing it with sequence lengths divisible by three
    * E.g. `atgaccatgattacgccaagcttgcatgcctgcagg[30,nt]` (encoding the 12-mer `MTMITPSLHACR` peptide) would produce oligos encoding every tiled 10-mer peptide.


* Ordinarily the script produces upper case output (e.g. `ttt	a[aa]	ccc` = `TTTGCCCCC`), but some options have been provided to more easily show the joins of the different input sequence regions.
  * The `-cf / --case_flip` flag makes sequential sequence fields alternate their case (e.g. `TTTgccCCC`)
  * The `--first_case` flag can be used in conjuction with this to specifically call the starting case (e.g. `-cf -fc lower` = `tttGCCccc`)


* There are also options to 3' pad nucleotides to ensure a minimum length (due to the requirements for certain oligo synthesis platforms)
  * This can be specified with the `-l / --oligo_len` flag, which will add random 3' nucleotides to the end of any oligos shorter than this.
  * Alternatively if further post-processing is to take place, additionally providing the `-n / --n_pad` flag will pad 3's with 'N' characters instead of random nucleotides.
