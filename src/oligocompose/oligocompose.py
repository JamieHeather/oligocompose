#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
oligocompose.py

Generate oligo pool library sequences from an input 'recipe' file
"""

import argparse
import collections as coll
import itertools as it
import pandas as pd
import gzip
import os
import random

__version__ = '0.2.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description="oligocompose v" + str(__version__) + '\n' +
                    ": Make large oligo pool sequences from an input 'recipe'")

    # Input and output options
    parser.add_argument('-in', '--in_file', required=True, type=str, default='',
                        help="Path to input recipe. Required.")

    parser.add_argument('-o', '--out_file', required=False, type=str, default='',
                        help="Path to output file. Optional. Default behaviour is 'out-[input-file-name].tsv'.")

    parser.add_argument('-f', '--fixed_file', required=False, type=str, default='',
                        help="Path to fixed sequence file. Optional.")

    parser.add_argument('-cf', '--case_flip', action='store_true', required=False, default=False,
                        help="Optional flag to make sequences from alternating columns flip case in the output.")

    parser.add_argument('-fc', '--first_case', type=str, required=False, default='upper',
                        help="Optional string ('upper' or 'lower') dictating the case of the first segment if "
                             "-cf/--case_flip flag invoked. Default = upper.")
    # TODO does this need a str parameter? -cf makes it upper case, maybe another flag to make it lower suffices?

    parser.add_argument('-l', '--oligo_len', type=int, required=False, default=0,
                        help="Integer of minimum length of oligos to produce (padding 3' with random nucleotides).")

    parser.add_argument('-n', '--n_pad', action='store_true', required=False, default=False,
                        help="Optional flag to overwrite random pad (-l) nucleotides with actual 'N' characters.")

    parser.add_argument('-nn', '--no_numbering', action='store_true', required=False, default=False,
                        help="Optional flag to prevent automatic numbering of output oligo sequences "
                             "(not recommended as it will likely result in many oligos with the same name).")

    parser.add_argument('-nl', '--number_length', type=int, required=False, default=4,
                        help="Integer of length of numbering to pad to (with leading zeroes). Default = 4.")

    parser.add_argument('-sl', '--step_length', type=int, required=False,
                        help="Step length for sliding windows. Defaults to 1 for AA (residues), 3 for NT (codons).")

    return parser.parse_args()


def get_seq_dict(path_to_file):
    """
    :param path_to_file: Path to file to be read in
    :return: Dict of named sequences
    """

    seq_dict = coll.defaultdict()
    with open(path_to_file, 'r') as in_file:
        for line in in_file:
            if line.startswith('#'):
                continue
            else:
                bits = fix(line).split('\t')
                if len(bits) == 2:
                    if bits[0] in seq_dict.keys():
                        raise IOError("Duplicate entry found in sequence file (" +
                                      path_to_file + "): " + bits[0])
                    else:
                        seq_dict[bits[0]] = bits[1].upper()
    return seq_dict


def opener(file_path, open_mode):
    """
    :param file_path: path to file to be opened
    :param open_mode: str detailing mode by which to open the file (e.g. w/r/a)
    :return: the appropriate file opening command (open or gzip.open)
    :raises: IOERrror: if file to be read doesn't exist, or if an unexpected file extension is specified
    """

    if open_mode == 'r':
        if not check_file_exists(file_path):
            raise IOError("Error, cannot locate file: " + file_path)

    if file_path.endswith('.gz'):
        return gzip.open(file_path, open_mode + 't')

    elif file_path.endswith('.tsv'):
        return open(file_path, open_mode)

    else:
        raise IOError("Error, unknown file extension: " + file_path)


def translate_nt(nt_seq):
    """
    :param nt_seq: Nucleotide sequence to be translated
    :return: corresponding amino acid sequence
    """

    aa_seq = ''
    for i in range(0, len(nt_seq), 3):
        codon = nt_seq[i:i+3].upper()
        if len(codon) == 3:
            try:
                aa_seq += codons[codon]
            except:
                raise IOError("Cannot translate codon: " + codon + ". ")

    return aa_seq


def rev_translate(amino_acid_seq, codon_usage_dict):
    """
    :param amino_acid_seq: An amino acid to convert into nucleotide sequence, using the most common codon
    :param codon_usage_dict: Dict of which codons to use for which amino acids (see get_optimal_codons)
    :return: Corresponding nucleotide sequence
    """

    return ''.join([codon_usage_dict[x] for x in amino_acid_seq])


def check_alphabet(possible_match, alphabet):
    """
    :param possible_match: A str of a sequence that may or may not derive from a particular alphabet (e.g. DNA/AA)
    :param alphabet: List of characters that constitute an 'alphabet'
    :return: True/False, relating to whether that sequence is plausible derived from that alphabet
    """
    return set(possible_match.upper()).issubset(list(alphabet))


def fix(string):
    """
    :param string: str to fix
    :return: sanitised str of 'string' with quote marks ('/") removed
    """
    return string.replace('"', '').replace("'", "").rstrip()


def check_file_exists(file_path):
    """
    :param file_path: str of path to a file
    :return: boolean as to whether that file exists
    """
    if os.path.exists(file_path):
        return True
    else:
        return False


def check_bioseq_type(bio_seq):
    """
    :param bio_seq: A str that may or may not be either a nucleotide or amino acid sequence
    :return: Str reflecting the inferred type
    """
    nt_chk = check_alphabet(bio_seq, ['A', 'C', 'G', 'T'])
    aa_chk = check_alphabet(bio_seq, opt_codons)
    if nt_chk:
        return 'dna'
    elif aa_chk:
        return 'protein'
    else:
        return 'X'


def get_seq_details(input_seq, key):
    """
    :param input_seq: str of a NT or AA sequence to process
    :param key: dict of parameter keys, allowing different functions to apply to different kinds of seqs
    :return: (str, dict), of updated output sequence and modified parameter key dict
    """

    start = input_seq.find('[')
    stop = input_seq.find(']')
    param_dict = coll.defaultdict()

    if start == -1 and stop == -1:
        out_seq = input_seq

    elif start != -1 and stop != -1:

        out_seq = input_seq[:start]
        fields = input_seq[start+1:stop].split(',')

        # First account for any potential numeric slice indications...
        param_dict['slice'] = []
        slices = [x for x in fields if x.isnumeric()]
        for seq_slice in slices:
            param_dict['slice'].append(int(seq_slice))
            fields.pop(fields.index(seq_slice))

        # ... including ranges
        slice_ranges = [x for x in fields if '-' in x]
        for srange in slice_ranges:
            sr_bits = srange.split('-')
            if len(sr_bits) != 2:
                raise ValueError("Improperly formatted slice range format (max one hyphen per field): " + srange)
            elif not [x.isnumeric() for x in sr_bits]:
                raise ValueError("Improperly formatted slice range format (only numbers allowed): " + srange)

            sr_bits = [int(x) for x in sr_bits]
            for sr in range(sr_bits[0], sr_bits[1] + 1):
                param_dict['slice'].append(sr)

            fields.pop(fields.index(srange))

        # Then check for the individual values
        for field in fields:
            if field.lower() in key.keys():
                param_dict[key[field]] = True
            else:
                raise ValueError("Unexpected sequence parameter detected: " + field)

    # Throw an error if a partial parameter tag detected
    elif start != -1 and stop == -1:
        raise ValueError("] not found in sequence " + input_seq)
    elif start == -1 and stop != -1:
        raise ValueError("] found without [ in sequence " + input_seq)

    # Determine if the 'sequence' is a reference
    if out_seq[0] in ['&', '$']:
        param_dict['reference'] = True
        if out_seq[0] == '&':
            param_dict['fixed'] = True
        elif out_seq[0] == '$':
            param_dict['file'] = True

    # Determine correct set of parameters given input
    if 'dna' in param_dict and 'protein' in param_dict:
        raise ValueError("Contradictory values provided ('dna' and 'protein') for sequence " + out_seq)

    elif 'dna' in param_dict and 'protein' not in param_dict:
        param_dict['type'] = 'dna'
    elif 'dna' not in param_dict and 'protein' in param_dict:
        param_dict['type'] = 'protein'

    else:
        inferred_type = check_bioseq_type(out_seq)
        if inferred_type != 'X':
            param_dict[inferred_type] = True
            param_dict['type'] = inferred_type
        elif 'reference' not in param_dict:
            raise ValueError("Unable to infer sequence type (nt/aa) for " + out_seq)

    if 'both' in param_dict:
        param_dict['forward'] = True
        param_dict['reverse'] = True
        param_dict.pop('both')

    if 'forward' not in param_dict and 'reverse' not in param_dict:
        param_dict['forward'] = True

    return out_seq, param_dict


def correct_case(string, desired_case):
    """
    :param string: A str relating to some sequence to ensure the correct case of
    :param desired_case: Str: 'upper' or 'lower'
    :return: Input string converted to the appropriate case
    """
    if desired_case.lower() == 'upper':
        return string.upper()
    elif desired_case.lower() == 'lower':
        return string.lower()
    else:
        raise IOError("Unrecognised case option detected: " + desired_case)


def compile_list_combinations(multi_size_list, final_output):
    """
    :param multi_size_list: list of sub-lists containing sequence combinations to be joined (e.g. [[a], [b,c], [d]])
    :param final_output: boolean of whether this is the final output (and thus whether case flipping can be attempted)
    :return: the conjoined combinations of those sub lists (e.g. [abd, acd])
    """

    # If case flip requested in input arguments, change the cases as appropriate
    if input_args['case_flip'] and final_output:
        for i in range(len(multi_size_list)):
            if i % 2 == 0:
                multi_size_list[i] = [correct_case(x, input_args['first_case']) for x in multi_size_list[i]]
            else:
                multi_size_list[i] = [correct_case(x, flip_case[input_args['first_case']]) for x in multi_size_list[i]]

    compiled = [''.join(x) for x in it.product(*multi_size_list)]

    return compiled


def expand_degeneracy(bio_seq, seq_type):
    """
    :param bio_seq: str of either DNA or amino acids
    :param seq_type: str detailing what the bio_seq is, either 'dna' or 'protein'
    :return: list of sequences
    """
    if seq_type == 'dna':
        conversion_dict = dna_codes
    elif seq_type == 'protein':
        conversion_dict = aa_codes
    else:
        raise ValueError("Invalid sequence type detected: " + seq_type)
    expanded_lists = [conversion_dict[x.upper()] for x in bio_seq]
    for x in compile_list_combinations(expanded_lists, False):
        yield x


def get_sliding_windows(seq_to_slide_across, window_len, step):
    """
    :param seq_to_slide_across: str of sequence from which to extract substrings
    :param window_len: int length of sliding window (and thus substrings produced)
    :param step: int step number in iteration, i.e. 1 for amino AA or 3 for DNA (assuming only coding codons wanted)
    :return: list of substrings of length window_len taken from seq_to_slide_across, at intervals of step
    """
    return [seq_to_slide_across[x:x + window_len]
            for x in range(0, len(seq_to_slide_across) - 1, step)
            if len(seq_to_slide_across[x:x + window_len]) == window_len]


def flatten_var_list(variable_content_list):
    out_list = []
    for item in variable_content_list:
        if isinstance(item, list):
            out_list.extend(flatten_var_list(item))
        else:
            out_list.append(item)
    return out_list


def parse_input_file(path_to_input_file, recursed, fixed):
    """
    :param path_to_input_file: str, path to file
    :param recursed: boolean, whether or not to allow recursion into files  # TODO keep ?
    :param fixed: boolean, whether or not there's a fixed sequence dict supplied
    :yield: processed sequences
    """

    if not check_file_exists(path_to_input_file):
        raise IOError("Can't find input file: " + path_to_input_file)

    with open(path_to_input_file, 'r') as in_file:

        line_count = 0
        for line in in_file:
            line_count += 1
            # Skip header/comment lines
            if line.startswith('#'):
                continue

            # Then determine what sequences need to be generated
            skip = False
            bits = fix(line).split('\t')
            nam = bits[0]
            combo_dict = coll.defaultdict(list)
            for section in range(1, len(bits)):

                seq, params = get_seq_details(bits[section], param_key)

                if 'reference' in params:

                    if 'fixed' in params:
                        seq_key = seq.replace('&', '')
                        if seq_key in fixed:
                            seq, new_params = get_seq_details(fixed[seq_key], param_key)
                            for np in new_params:
                                if np not in params:
                                    params[np] = new_params[np]
                        else:
                            raise ValueError("Specific reference sequence not found: " + seq)

                    # TODO leave this for now, think it might have to be a separate function
                    if 'file' in params:
                        if recursed:
                            print("Warning: detected a reference to another file (" + seq + "): other files may "
                                  "only be referenced from the top-level oligo tsv file - skipping.")
                            skip = True
                            continue
                        path_to_file = seq.replace('$', '')
                        if check_file_exists(path_to_file):
                            seq = parse_input_file(path_to_file, True, fixed)

                # If the sequence(s) need slicing, chop 'em up
                if 'slice' in params:
                    if params['slice']:
                        if input_args['step_length']:
                            window_step = input_args['step_length']
                        elif params['type'] == 'dna':
                            window_step = 3
                        elif params['type'] == 'protein':
                            window_step = 1

                        seqs = []
                        for s in params['slice']:
                            seqs.append(get_sliding_windows(seq, s, window_step))
                        seqs = flatten_var_list(seqs)
                    else:
                        seqs = [seq]
                else:
                    seqs = [seq]

                seqs = flatten_var_list(seqs)

                # Expand degenerate bases
                for seq_i in range(len(seqs)):
                    seqs[seq_i] = [x for x in expand_degeneracy(seqs[seq_i], params['type'])]

                seqs = flatten_var_list(seqs)

                # Reverse translate those as need it
                if params['type'] == 'protein':
                    for seq_i in range(len(seqs)):
                        seqs[seq_i] = rev_translate(seqs[seq_i], opt_codons)

                if skip:
                    continue

                # Finally stick all the bits together
                combo_dict[section] = (seqs, params)

            compiled_seqs = compile_list_combinations([combo_dict[x][0] for x in combo_dict], True)

            if input_args['no_numbering']:
                yield [(nam, compiled_seqs[x]) for x in range(len(compiled_seqs))]

            else:
                yield [(nam + '|' + str(x + 1).zfill(input_args['number_length']),
                        compiled_seqs[x]) for x in range(len(compiled_seqs))]


def list_to_df(input_list, headers, rename):
    """
    Convert a list to a (long) dataframe. Note that first entry becomes the index if chosen
    :param input_list: List of list entries (with each position in each list corresponding to a column)
    :param headers: List of column headers. First column should be unique, becoming the rownames, if rename = True
    :param rename: Option to rename row IDs by first colum
    :return: sorted pandas dataframe
    """
    out_df = pd.DataFrame(input_list)
    out_df = out_df.rename(index=str, columns=dict(zip(range(len(headers)), headers)))
    out_df = out_df.sort_values(by=[headers[0]])
    if rename is True:
        out_df = out_df.set_index(headers[0], drop=True)
    return out_df


def random_pad(length_to_pad, avoid_list, case='lower'):
    """
    :param length_to_pad: int, number of random nucleotides to produce
    :param avoid_list: list of str, sequences that aren't allowed in the final random sequence
    :param case: str, case of output sequence, 'lower' = default, anything else produces upper case
    :return: a random DNA string of a defined length (with every base being equally likely)
    """
    avoid_list += ['aaaa', 'cccc', 'gggg', 'tttt', 'AAAA', 'CCCC', 'GGGG', 'TTTT']
    clash_free = False
    while not clash_free:
        pad = ''.join([random.choice('ACGT') for x in range(length_to_pad)])
        clash_check = [x for x in avoid_list if x in pad]
        if len(clash_check) == 0:
            clash_free = True
    if case == 'lower':
        return pad.lower()
    else:
        return pad.upper()


param_key = {
    'f': 'forward',
    'r': 'reverse',
    'b': 'both',
    'nt': 'dna',
    'aa': 'protein',
    'o': 'optimised',
    'p': 'provided',
    's': 'scan'         # TODO provide option to do alanine-/x-scans?
}

dna_codes = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], 'U': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['C', 'G'],
    'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

all_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

flip_case = {'upper': 'lower', 'lower': 'upper'}

aa_codes = {
    'A': ['A'], 'C': ['C'], 'D': ['D'], 'E': ['E'],
    'F': ['F'], 'G': ['G'], 'H': ['H'], 'I': ['I'],
    'K': ['K'], 'L': ['L'], 'M': ['M'], 'N': ['N'],
    'P': ['P'], 'Q': ['Q'], 'R': ['R'], 'S': ['S'],
    'T': ['T'], 'V': ['V'], 'W': ['W'], 'Y': ['Y'],
    'X': all_aa, '!': all_aa + ['*']
}

codons = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
          'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
          'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
          'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
          'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
          'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
          'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
          'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
          'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
          'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
          'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
          'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
          'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
          'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
          'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
          'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}

# Most commonly used human codons (from Kazusa via Stitchr)
opt_codons = {'F': 'TTC', 'S': 'AGC', 'Y': 'TAC', 'C': 'TGC', 'L': 'CTG', '*': 'TGA', 'W': 'TGG',
              'P': 'CCC', 'H': 'CAC', 'R': 'AGA', 'Q': 'CAG', 'I': 'ATC', 'T': 'ACC', 'N': 'AAC',
              'K': 'AAG', 'M': 'ATG', 'V': 'GTG', 'A': 'GCC', 'D': 'GAC', 'G': 'GGC', 'E': 'GAG'}


def main():

    global input_args
    input_args = vars(args())

    # Read in fixed sequences
    if input_args['fixed_file']:
        fixed_dict = get_seq_dict(input_args['fixed_file'])
    else:
        fixed_dict = ''

    # Determine output filename
    if input_args['out_file']:
        if check_file_exists(input_args['out_file']) or check_file_exists(input_args['out_file'] + '.tsv'):
            print("Overwriting existing outfile")
    else:
        input_args['out_file'] = 'out-' + input_args['in_file'].split('/')[-1]

    # Then go through the requested file and generate the sequences
    data = flatten_var_list(parse_input_file(input_args['in_file'], False, fixed_dict))

    # Convert that to a dataframe, merge redundant rows
    df = list_to_df(data, ['ID', 'Sequence'], False)
    df = df.groupby('Sequence', as_index=False).aggregate({'ID': '_'.join}).reindex(columns=df.columns)
    df = df.sort_values('ID')

    # If padding needed, get padding
    if input_args['oligo_len'] > 0:
        padded_seqs = []
        for i in range(len(df)):
            seq = df.iloc[i]['Sequence']
            pad_difference = input_args['oligo_len'] - len(seq)

            if pad_difference > 0:

                if input_args['n_pad']:
                    if input_args['case_flip']:
                        if seq[-1].isupper():
                            padded_seqs.append(seq + 'n' * pad_difference)
                        else:
                            padded_seqs.append(seq + 'N' * pad_difference)

                else:
                    pad_case = 'upper'  # Default
                    if input_args['case_flip']:
                        if seq[-1] == seq[-1].lower():
                            pad_case = 'upper'
                        else:
                            pad_case = 'lower'

                    # TODO engineer in avoid list option
                    padded_seqs.append(seq + random_pad(pad_difference, [], pad_case))
        df['Sequence'] = padded_seqs

    # Finally output
    df.to_csv(input_args['out_file'], sep='\t', index=False, header=False)
    print("Output " + str(len(df)) + " oligos to " + input_args['out_file'])
