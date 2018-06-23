from collections import defaultdict, Counter, Iterable
from itertools import islice, product, chain
from warnings import warn

import numpy as np
import Bio.Data.CodonTable
from Bio.Seq import Seq

def create_id_from_row(row, cols_to_concatenate):

    id = ''

    for c in cols_to_concatenate:

        z = row[c]
        z = z.replace('"', '')
        z = z.replace(';', '')

        id += z + '_'

    return id

# want to kmerize long rna sequences
# make a function that takes in a subsequence and a k, returns a vector of kmer frequences

def k_mers(seq, k):
    '''Yields all *k*-mers in the input sequence with repeats.
    Args:
        seq (str): The sequence for which to generate *k*-mers.
        k (int): the length of the *k*-mers.
    Yields:
        str: the next *k*-mer
    Raises:
        ValueError: When the value of *k* is less than the length of the sequence, k <= 0, or len(seq) is 0.
    Example:
        >>> list(k_mers("GATTACA", 1))
        ['G', 'A', 'T', 'T', 'A', 'C', 'A']
        >>> list(k_mers("GATTACA", 2))
        ['GA', 'AT', 'TT', 'TA', 'AC', 'CA']
        >>> list(k_mers("GATTACA", 3))
        ['GAT', 'ATT', 'TTA', 'TAC', 'ACA']
        >>> list(k_mers("GATTACA", 4))
        ['GATT', 'ATTA', 'TTAC', 'TACA']
        >>> k_mers("GATTACA", 4)
        <generator object k_mers at 0x10831d258>
    '''

    # error checking
    if k > len(seq):
        raise ValueError("k (%i) may not be less then length of seq (%i)." % (k, len(seq)))
    elif len(seq) == 0:
        raise ValueError("seq length may not be zero")
    elif k <= 0:
        raise ValueError("k may not be <= zero")

    it = iter(seq)
    result = tuple(islice(it, k))
    if len(result) == k:
        yield "".join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield "".join(result)

def k_mer_frequencies(seq, k, include_missing=True, vector=False):
    '''Calculates relative frequencies of each *k*-mer in the sequence.
    Args:
        seq (str or list): The sequence(s) to for which to generate *k*-mer frequencies.
        k (int or list): the length of the *k*-mer(s).
        include_missing (bool, optional): If True, include missing *k*-mers as having a frequency of 0. Only supports DNA *k*-mers. Defaults to False.
        vector (bool, optional): Return a 1-D Numpy array of the *k*-mer frequencies, ordered by *k*-mers alphabetically. If True, ``include_missing`` must also be True. Defaults to False.
    Returns:
        dict: A dict in which the keys are *k*-mers and the values are floats of their frequencies.
    Raises:
        ValueError: When an invalid value of k is provided or ``include_missing`` is False and ``vector`` is True.
        ValueError: When ``k`` or ``seq`` is not provided.
    Example:
        >>> k_mer_frequencies("INQTEL", 1, include_missing=False)
        {'E': 0.16666666666666666,
         'I': 0.16666666666666666,
         'L': 0.16666666666666666,
         'N': 0.16666666666666666,
         'Q': 0.16666666666666666,
         'T': 0.16666666666666666}
        >>> k_mer_frequencies("GATGATGGC", [1, 2], include_missing=False)
        {'A': 0.2222222222222222,
         'AT': 0.25,
         'C': 0.1111111111111111,
         'G': 0.4444444444444444,
         'GA': 0.25,
         'GC': 0.125,
         'GG': 0.125,
         'T': 0.2222222222222222,
         'TG': 0.25}
        >>> k_mer_frequencies(["A", "T"], 1, include_missing=False)
        {"A": 0.5, "T": 0.5}
        >>> k_mer_frequencies("GATGATGGC", 2, include_missing=True)
        {'AA': 0,
         'AC': 0,
         'AG': 0,
         'AT': 0.25,
         'CA': 0,
         'CC': 0,
         'CG': 0,
         'CT': 0,
         'GA': 0.25,
         'GC': 0.125,
         'GG': 0.125,
         'GT': 0,
         'TA': 0,
         'TC': 0,
         'TG': 0.25,
         'TT': 0}
        >>> k_mer_frequencies("GATGATGGC", 2, include_missing=True, vector=True)
        array([0.   , 0.   , 0.   , 0.25 , 0.   , 0.   , 0.   , 0.   , 0.25 ,
               0.125, 0.125, 0.   , 0.   , 0.   , 0.25 , 0.   ])
    '''

    if not include_missing and vector:
        raise ValueError("May not create vector without including missing kmers.")
    elif not k:
        raise ValueError("Must provide a value for k")
    elif not seq:
        raise ValueError("Must provide seq(s)")

    if not isinstance(k, Iterable):
        k = [k]
    else:
        k = sorted(k)

    output = []

    if isinstance(seq, (str, bytes, Seq)):
        seq = [seq]

    for _k in k:

        # check the value of k
        if _k < 1:
            raise ValueError("Invalid value of k. May not be less than 1.")

        # get all the k-mers for the seqs
        _seqs = []
        for _seq in [list(k_mers(_seq, _k)) for _seq in seq]:
            _seqs.extend(_seq)

        # determine their frequencies
        count = Counter(_seqs)
        frequencies = {key: value / sum(count.values()) for key, value in count.items()}

        if include_missing:
            defaults = {"".join(x): 0 for x in list(product("ATGC", repeat=_k))}
            frequencies = {**defaults, **frequencies}
        if vector:
            frequencies = sorted(list(frequencies.items()), key=lambda x: x[0])
            frequencies = np.fromiter((x[1] for x in frequencies), float, count=len(frequencies))
        output.append(frequencies)

    if len(output) == 1:
        return output[0]
    elif vector:
        return np.array(list(chain.from_iterable(output)))
    else:
        return {k: v for d in output for k, v in d.items()}

#print(k_mer_frequencies('ATCGGGATCTGCAGAGTTTCAGTCGATCGATCGGCTAGCAAACTA', [1,2,4], True, True))

# run this on each row, generate a vector of k_mer_frequencies
# train a model on these feature vectors with the overlap labels
# divide each tissue type into train and test
