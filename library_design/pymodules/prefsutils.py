"""
This script contains functions for analylzing amino acid preferences

Hugh Haddox, February-19-2016
"""

import re
import math
import doctest
import subprocess
import pandas
import scipy.stats
import scipy.optimize
import numpy

def amino_acids(stop=False):
    """This function returns the amino acids. If stop==True, it returns a stop character '*' as the 21st item"""
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    assert len(amino_acids) == 20
    if stop:
        amino_acids += ['*']
        assert len(amino_acids) == 21
    return amino_acids

def list_all_single_nt_changes_to_codon(WTcodon):
    """
    This function returns all single-nucleotide changes to a given wildtype codon
    
    >>> list_all_single_nt_changes_to_codon('ATA')
    ['AAA', 'TTA', 'ATT', 'CTA', 'ACA', 'ATC', 'GTA', 'AGA', 'ATG']
    """
    WTcodon = WTcodon.upper()
    assert len(WTcodon) == 3
    
    SingleNtChanges = []
    for nt in ['A', 'T', 'C', 'G']:
        CodonVariantPosition0 = nt + WTcodon[1] + WTcodon[2]
        CodonVariantPosition1 = WTcodon[0] + nt + WTcodon[2]
        CodonVariantPosition2 = WTcodon[0] + WTcodon[1] + nt
        if CodonVariantPosition0 != WTcodon:
            SingleNtChanges.append(CodonVariantPosition0)
        if CodonVariantPosition1 != WTcodon:
            SingleNtChanges.append(CodonVariantPosition1)
        if CodonVariantPosition2 != WTcodon:
            SingleNtChanges.append(CodonVariantPosition2)
    
    assert len(SingleNtChanges) == 9
    
    return SingleNtChanges

def ComputeEntropy(frequencies, base):
    """
    This function computes the entropy of a list of frequencies (*frequencies*) for a given
    log base (*base*)
    
    Code for doctest:
    >>> base = 2
    
    >>> frequencies = [1.0, 1e-11, 0.0]
    >>> ComputeEntropy(frequencies, base)
    0.0
    
    >>> frequencies = [0.5, 0.5]
    >>> ComputeEntropy(frequencies, base)
    1.0
    
    >>> frequencies = [0.2, 0.8]
    >>> '%.6f'%ComputeEntropy(frequencies, base)
    '0.721928'
    """
    
    if abs(1 - sum(frequencies)) > 1e-5:
        raise ValueError("The frequencies do not add up to 1 +/- 1e-5")
    
    entropy = sum([(-f*math.log(f, base)) for f in frequencies if f > 1e-10])
    
    return entropy

def ComputeSimpsonsDiversity(frequencies):
    """
    This function computes the Simpson's diversity of a list of frequencies (*frequencies*)
    
    Code for doctest:
    
    >>> frequencies = [1.0, 1e-11, 0.0]
    >>> ComputeSimpsonsDiversity(frequencies)
    1.0
    
    >>> frequencies = [0.5, 0.5]
    >>> ComputeSimpsonsDiversity(frequencies)
    0.5
    
    >>> frequencies = [0.2, 0.8]
    >>> '%.2f'%ComputeSimpsonsDiversity(frequencies)
    '0.68'
    """
    diversity = sum([math.pow(f, 2) for f in frequencies])
    return diversity

def TranslateCodon(codon):
    """Returns one-letter amino acid code for *codon*.

    *codon* is a 3-letter string giving a valid codon."""
    genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
        'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
        'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
        'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
        'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
        'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
        'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
        'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
        'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W', 'CGT':'R',
        'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
        'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    return genetic_code[codon.upper()]

def write_sbatch_file(cmd, sbatchfilename, ncpus, mem=50000):
    """This function writes an sbatch file using the following input variables:
            *cmd*: is the command to be executed
            *sbatchfilename*: is the name of the sbatch file and the prefix for outfiles/error files"""
    with open(sbatchfilename, 'w') as sbatchfile:
        sbatchfile.write("#!/bin/bash\n")
        sbatchfile.write("#SBATCH\n")
        sbatchfile.write("#SBATCH -o {0}.out\n".format(sbatchfilename))
        sbatchfile.write("#SBATCH -e {0}.err\n".format(sbatchfilename))
        sbatchfile.write("#SBATCH -p largenode\n")
        sbatchfile.write("#SBATCH -c {0}\n".format(ncpus))
        sbatchfile.write("#SBATCH -t 1-00:00\n")
        sbatchfile.write("#SBATCH --mem={0}\n".format(mem))
        sbatchfile.write(cmd)
    pass

def ShannonJensenDivergence(p1, p2):
    """Returns the Shannon-Jensen divergence between *p1* and *p2*.

    Requires ``numpy``, which is available if *PymcAvailable()* evaluates
    to *True*. Otherwise you will get an error.

    *p1* and *p2* are two *numpy.ndarray* objects with elements that give
    the probability distributions for which we are computing the divergence.
    Specifically, *p1[i]* and *p2[i]* indicate the probability for state *i*
    in the two distributions. The constraints are:
    
        * *p1* and *p2* should have dimension of one
          (*1 == len(p1.shape) == len(p2.shape)*).
          
        * *p1* and *p2* should have the same
          nonzero length (*len(p1) == len(p2) > 0*)
          
        * *p1* and *p2* should both have entries that sum to one (since
          they are probability distributions), so
          *True == numpy.allclose(sum(p1), 1) == numpy.allclose(sum(p2), 1)*

        * *p1* and *p2* should have all entries >= 0, so 
          *True == numpy.all(p1 >= 0) == numpy.all(p2 >= 0)*.

    This function returns the Shannon-Jensen divergence between the 
    probability distributions specified by *p1* and *p2*. The logarithms
    are taken to the base 2, meaning that the returned divergence
    will always be between 0 and 1.
    """
    assert 1 == len(p1.shape) == len(p2.shape), "p1 and p2 not both numpy.ndarrays of dimension 1"
    assert len(p1) == len(p2) > 0, "p1 and p2 not both arrays of same nonzero length"
    assert numpy.allclose(sum(p1), 1, atol=0.005), "p1 does not have entries summing to one: sum is %g, p1 is %s" % (sum(p1), p1)
    assert numpy.allclose(sum(p2), 1, atol=0.005), "p2 does not have entries summing to one: sum is %g, p2 is %s" % (sum(p2), p2)
    assert numpy.all(p1 >= 0), "p1 does not have all entries >= 0: p1 is %s" % p1
    assert numpy.all(p2 >= 0), "p2 does not have all entries >= 0: p2 is %s" % p2
    m = (p1 + p2) / 2.0
    d_p1_m = d_p2_m = 0.0
    for i in range(len(p1)):
        if p1[i]:
            d_p1_m += p1[i] * math.log(p1[i] / m[i], 2)
        if p2[i]:
            d_p2_m += p2[i] * math.log(p2[i] / m[i], 2)
    jsd = (d_p1_m + d_p2_m) / 2.0
    assert -1e10 <= jsd <= 1, "Shannon-Jensen divergence should be between zero and one, got value of %g" % jsd
    return jsd

def ComputeDistanceBetweenPreferencesAtASite(prefs1, prefs2, site, distance_metric, stop=True):
    """
    This function computes the distance between the preferences of two homologs at a given site. Given the following input:
    *prefs1* and *prefs2* : dictionaries of site-specific amino-acid preferences for homologs 1 and 2
    *site* : site for consideration
    *distance_metric* : distance metric to be used. See below for options:
            * *half_sum_of_absolute_diffs* : subtracts two vectors of preferences, sums the absolute value of each value
            in the resulting vector, and then divides that value by two so that the range of the final value is between 0-1
            * *normalized_RMSD* : computes the RMSD between two vectors of preferences and multiplies this by the square root of
            1/2 to normalize the values such that they range between 0-1.
            * *Jensen_Shannon_distance* : the square root of the Jensen-Shannon divergence.   
            * *hydropathy* : the absolute value of the difference in site-specific hydropathy scores
    *stop* : a boolean specifying whether to include the preferences for stop codons in the analysis
    
    The output of this function is the computed distance.
    
    Code for doctest:
    >>> site = '165'
    >>> prefs1 = {'165': {'WT': 'L', 'entropy': 2.4851505880264964, '*': 0.00404600492090634, 'A': 0.012360535030154038, 'C': 0.04229213916539216, 'E': 0.01735503332858423, 'D': 0.016008205436725518, 'G': 0.009273300423732627, 'F': 0.013103952684284674, 'I': 0.03416518082804878, 'H': 0.010859315515897245, 'K': 0.012875958082953168, 'M': 0.018343904825321222, 'L': 0.6253814390341444, 'N': 0.005943697926058552, 'Q': 0.011095374141738118, 'P': 0.011522399843107588, 'S': 0.01702395912623408, 'R': 0.019388049942067616, 'T': 0.024661437843780845, 'W': 0.023837772699945227, 'V': 0.021228303965178216, 'Y': 0.0492340352357455}}
    >>> prefs2 = {'165': {'WT': 'I', 'entropy': 3.8413342561337807, '*': 0.01977055310212826, 'A': 0.06667787718773868, 'C': 0.0199059706587126, 'E': 0.048803607638508904, 'D': 0.0475664057425134, 'G': 0.008085134144643064, 'F': 0.0049025553364354605, 'I': 0.09106641456379641, 'H': 0.012503576115388321, 'K': 0.05186996939263794, 'M': 0.11069490149040916, 'L': 0.09949104954328204, 'N': 0.009046382157482687, 'Q': 0.029057404698214426, 'P': 0.10011977750934838, 'S': 0.019640847549499355, 'R': 0.01635178548146221, 'T': 0.0023214534983395188, 'W': 0.037801646527372036, 'V': 0.030955756082615336, 'Y': 0.17336693157947153}}
    >>> round(ComputeDistanceBetweenPreferencesAtASite(prefs1, prefs2, site, 'half_sum_of_absolute_diffs'), 5)
    0.58304
    >>> round(ComputeDistanceBetweenPreferencesAtASite(prefs1, prefs2, site, 'normalized_RMSD'), 6)
    0.399999
    """
    if stop:
        aas = amino_acids(stop=True)
    else:
        aas = amino_acids(stop=False)
    site_prefs_1 = [prefs1[site][aa] for aa in aas]
    site_prefs_2 = [prefs2[site][aa] for aa in aas]
    if distance_metric == 'half_sum_of_absolute_diffs':
        distance = sum([abs(p1-p2) for (p1,p2) in zip(site_prefs_1, site_prefs_2)]) / 2.0
    elif distance_metric == 'normalized_RMSD':
        distance = math.sqrt( sum([math.pow(p1-p2,2) for (p1,p2) in zip(site_prefs_1, site_prefs_2)]) / 2.0 )
    elif distance_metric == 'Jensen_Shannon_distance':
        distance = math.sqrt(ShannonJensenDivergence(numpy.asarray(site_prefs_1), numpy.asarray(site_prefs_2)))
    elif distance_metric == 'hydropathy':
        site_hydropathies_1 = CalculateAAPreferenceHydropathy(prefs1)
        site_hydropathies_2 = CalculateAAPreferenceHydropathy(prefs2)
        distance = abs(site_hydropathies_1[site] - site_hydropathies_2[site])
    else:
        raise ValueError('Did not recognize distance metric: {0}'.format(distance_metric))
    return distance
        
def ComputeRMSD(vector):
    """
    This function computes the root mean square distance of a vector of values (*vector*).
    Code for doctest:
    >>> v = [1.2, 3.5, 6.8, 1.1]
    >>> round(ComputeRMSD(v), 5)
    3.9096
    """
    RMSD = math.sqrt( sum([math.pow(float(d),2) for d in vector]) / float(len(vector)) )
    return RMSD

if __name__ == '__main__':
    import doctest
    doctest.testmod()