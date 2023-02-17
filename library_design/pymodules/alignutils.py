"""
This script contains functions for working with FASTA files.

Note: many of these functions assume alignments are in the form of a list of tupples in the form:
    [(header1, seq1), (header2, seq2), ...]

Hugh Haddox, February-19-2016
"""

from Bio import SeqIO
import re
import math
import doctest
import subprocess


def ReadInAlignment(alignment_file):
    """This function reads in an alignment in FASTA form from a given file (*alignment_file*). It then
    returns the alignment in the form of a list: [(header1, seq1), (header2, seq2), ...]"""
    a = []
    with open(alignment_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            (header, seq) = (record.id, str(record.seq))
            a.append((header, seq))
    return a


def WriteAlignmentToFASTA(a, outfile):
    """This function writes an alignment *a* to a file *outfile* in FASTA format"""
    with open(outfile, 'w') as f:
        for (header, seq) in a:
            f.write('>{0}\n{1}\n'.format(header, seq))
    return None


def AssertHeadersMatch(headers_seqs_1, headers_seqs_2):
    """This function asserts that two alignments (*headers_seqs_1* and *headers_seqs_2* in the list form
    described above) have the exact same headers in the same order"""
    assert len(headers_seqs_1) == len(headers_seqs_2), "The alignments have different lengths: {0} and {1}".format(len(headers_seqs_1), len(headers_seqs_2))
    for i in range(len(headers_seqs_1)):
        if headers_seqs_1[i][0] != headers_seqs_2[i][0]:
            raise ValueError("The headers for entry %s don't match: %s %s"%(i+1, headers_seqs_1[i][0], headers_seqs_2[i][0]))
    pass


def FilterAlignmentbyLengthandDisallowedCharacters(headers_seqs, disallowed_characters, length):
    """
    This function removes sequences from an alignment (*headers_seqs*; in the list form described above)
    if they don't match the specified *length* (int), if they contain any of the specified
    *disallowed_characters* (list), or if they have premature stop codons.
    
    Code for doctest:
    >>> disallowed_characters = ['#', 'X']
    >>> test_headers_seqs = [('s1', 'A-G-T*'), ('s2', 'AGG-A'), ('s3', 'AGG-A*'), ('s4', 'A-#-T*'), ('s5', 'X-G-T*'), ('s6', 'A*G-T*')] 
    >>> length = len(test_headers_seqs[0][1])
    >>> (a, d) = FilterAlignmentbyLengthandDisallowedCharacters(test_headers_seqs, disallowed_characters, length)
    >>> a
    [('s1', 'A-G-T*'), ('s3', 'AGG-A*')]
    >>> d
    [('s2', 'AGG-A'), ('s4', 'A-#-T*'), ('s5', 'X-G-T*'), ('s6', 'A*G-T*')]
    """
    discarded_seqs = []
    filtered_headers_seqs = []
    for (header, seq) in headers_seqs:
        discard_seq = False
        seq = seq.upper()
        # Filter sequences with premature stop codons
        if '*' in seq[:-1]:
            discard_seq = True
        # Filter sequences with disallowed characters
        for char in disallowed_characters:
            if char in seq:
                discard_seq = True
        # Filter sequences that don't match the specified length
        if len(seq) != length:
            discard_seq = True
        # Append to the filtered alignment sequences that passed the filtering steps
        if discard_seq == False:
            filtered_headers_seqs.append((header, seq))
        else:
            discarded_seqs.append((header,seq))
    assert len(headers_seqs) == len(filtered_headers_seqs)+len(discarded_seqs)
    return (filtered_headers_seqs, discarded_seqs)


def StripGapsToFirstSequence(a):
    """This function takes an alignment *a* and removes all columns where the first sequence has a gap.
    Test code for doctest:
    >>> a = [('seq1', 'AT-CG-A'), ('seq2', 'ATG-GGA')]
    >>> StripGapsToFirstSequence(a)
    [('seq1', 'ATCGA'), ('seq2', 'AT-GA')]
    """
    seq1 = a[0][1]
    seq1_gaps = [i for (i,s) in enumerate(seq1) if s=='-']
    new_a = []
    for (header, seq) in a:
        new_a.append((header, ''.join([s for (i,s) in enumerate(seq) if i not in seq1_gaps])))
    return new_a


def IdentifyColumnsInAlignmentWithNumberOfGapsAboveThreshold(headers_seqs, f_threshold):
    """
    This function identifies columns in an alignment where the number of gaps in
    that column is greater than some threshold.
    
    Code for doctest:
    >>> test_headers_seqs = []
    >>> test_headers_seqs.append(('s1', 'A-G-T*'))
    >>> test_headers_seqs.append(('s2', 'AGG--*'))
    >>> test_headers_seqs.append(('s3', 'AC---*'))
    >>> test_headers_seqs.append(('s4', 'AG---*'))
    >>> IdentifyColumnsInAlignmentWithNumberOfGapsAboveThreshold(test_headers_seqs, 0.51)
    [4, 5]
    >>> IdentifyColumnsInAlignmentWithNumberOfGapsAboveThreshold(test_headers_seqs, 0.2)
    [2, 3, 4, 5]
    """
    # I will count the number of gaps in each column, with columns indexed starting at 1
    seq_len = len(headers_seqs[0][1])
    ngaps = dict((i+1, 0) for i in range(seq_len))
    for (header, seq) in headers_seqs:
        for i in ngaps:
            if seq[i-1] == '-':
                ngaps[i] += 1
    # Next, I will make a list of columns for which the number of gaps exceeds the threshold
    n_threshold = float(f_threshold) * float(len(headers_seqs))
    columns_above_threshold = []
    for i in ngaps:
        if ngaps[i] > n_threshold:
            columns_above_threshold.append(i)
    columns_above_threshold.sort()
    return columns_above_threshold


def CreateDNAColumnListFromProteinColumnList(protein_columns):
    """This function takes a list of columns in a protein alignment and returns the corresponding list
    of columns for a DNA alignment of the same gene
    
    Code for doctest:
    >>> test_protein_columns = [1, 3, 8]
    >>> CreateDNAColumnListFromProteinColumnList(test_protein_columns)
    [1, 2, 3, 7, 8, 9, 22, 23, 24]
    """
    dna_columns = []
    for protein_column in protein_columns:
        first_dna_column_in_codon = (3*protein_column)-2
        dna_columns.extend([first_dna_column_in_codon, first_dna_column_in_codon + 1, first_dna_column_in_codon + 2])
    return dna_columns


def FilterAlignmentByPatternInHeader(headers_seqs, re_pattern):
    """
    This function takes sequences from a FASTA file and only retains those that have the specified pattern
    (string that will be converted into a regular expression) in their header.
    
    Code for doctest:
    >>> test_headers_seqs = [('a.5.s1', 'A-G-T*'), ('s2', 'AGG--*'), ('s3a.5', 'ACTT-*'), ('a.5.s4', 'AG---*')]
    >>> test_pattern = '^a\.5'
    >>> FilterAlignmentByPatternInHeader(test_headers_seqs, test_pattern)
    [('a.5.s1', 'A-G-T*'), ('a.5.s4', 'AG---*')]
    """
    pattern = re.compile('%s'%re_pattern)
    retained_headers_seqs = [(header, seq) for (header, seq) in headers_seqs if re.search(pattern, header)]
    return retained_headers_seqs


def RemoveColumnsFromAlignment(headers_seqs, columns_for_removal):
    """
    This function removes columns from an alignment
    
    Code for doctest:
    >>> test_headers_seqs = [('s1', 'A-G-T*'), ('s2', 'AGG--*'), ('s3', 'ACTT-*'), ('s4', 'AG---*')]
    >>> RemoveColumnsFromAlignment(test_headers_seqs, [1, 4, 6])
    [('s1', '-GT'), ('s2', 'GG-'), ('s3', 'CT-'), ('s4', 'G--')]
    """
    shortened_headers_seqs = []
    for (header, seq) in headers_seqs:
        shortened_seq = []
        for i in range(len(seq)):
            if i+1 not in columns_for_removal:
                shortened_seq.append(seq[i])
        assert len(seq)-len(columns_for_removal) == len(shortened_seq)
        shortened_headers_seqs.append((header, ''.join(shortened_seq)))
    
    return (shortened_headers_seqs)


def Create_Remove_And_Renumber_Files_For_dms_editsites(columns_for_removal, removefile, renumberfile, HXB2_bounds):
    """
    This function creates two files for *dms_editsites*. These files are for removing and renumbering
    sites in .txt files of preferences so that they match alignments after removal of a subset of columns:
    *columns_for_removal* : a list of strings specifying which columns should be removed. Note: this should include
    non-HXB2 sites for a given homolog.
    *removefile* : a file that specifies sites for removal.
    *renumberfile* : a file that specifies sites for renumbering.
    *HXB2_bounds* : a tupple giving the bounds of the preferences files in HXB2 numbering (integers; inclusive)
    """
    # First, I will make a file specifying which sites in HXB2 to remove.
    with open(removefile, 'w') as dms_remove:
        dms_remove.write('# sites to remove\n')
        for site in columns_for_removal:
            dms_remove.write('%s\n'%site)

    # Next, I will make a file to renumber the remaing sites in sequential order using 1, 2, ... numbering.
    with open(renumberfile, 'w') as dms_renumber:
        dms_renumber.write('#ORIGINAL_SITE NEW_SITE\n')
        n = 1 # indexed starting at 1
        (low_bound, high_bound) = HXB2_bounds
        for site in range(low_bound, high_bound+1):
            if str(site) not in columns_for_removal:
                dms_renumber.write('%s %s\n'%(site, n))
                n += 1

def ComputePercentSequenceIdentity(seq1, seq2):
    """
    This function computes the percent sequence identity for two sequences (strings; case insensitive). Positions
    where both sequences have a gap character are ignored.
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    assert len(seq1) == len(seq2)
    matches_bools = [i==j for (i,j) in zip(seq1, seq2) if not i==j=='-']
    n_sites_compared = len(matches_bools)
    percent_ID = 100 * float(sum(matches_bools))/n_sites_compared
    return (percent_ID, n_sites_compared)
                
if __name__ == '__main__':
    import doctest
    doctest.testmod()