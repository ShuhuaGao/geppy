"""
.. module:: geppy.tools.mutation
.. moduleauthor:: Shuhua Gao

The module :mod:`mutation` provides mutation related genetic modifications in GEP, including point-mutation,
transposition and inversion. Please refer to Chapter 3 of [FC2006]_ for more details.

"""
import random
from ..core.symbol import Function

_DEBUG = False


def _choose_function(pset):
    return random.choice(pset.functions)


def _choose_terminal(pset):
    return random.choice(pset.terminals)


def _choose_subsequence(seq, min_length=1, max_length=-1):
    if max_length <= 0:
        max_length = len(seq)
    length = random.randint(min_length, max_length)
    start = random.randint(0, len(seq) - length)
    return start, start + length


def mutUniform(individual, pset, indpb):
    """
    For each symbol in *individual*, change it to another randomly chosen symbol from *pset* with the probability
    *indpb*.

    :param individual: the chromosome to be mutated.
    :param pset: :class:`PrimitiveSet`, a primitive set
    :param indpb: probability of mutating each symbol
    :return: A tuple of one chromosome

    It is typical to set a mutation rate *indpb* equivalent to two one-point mutations per chromosome.
    """
    for gene in individual:
        # mutate the gene with the associated pset
        # head: any symbol can be changed into a function or a terminal
        for i in range(gene.head_length):
            if random.random() < indpb:
                if random.random() < 0.5:  # to a function
                    gene[i] = _choose_function(pset)
                else:
                    gene[i] = _choose_terminal(pset)
        # tail: only change to another terminal
        for i in range(gene.head_length, len(gene)):
            if random.random() < indpb:
                gene[i] = _choose_terminal(pset)
    return individual,


def invert(individual):
    """
    A gene is randomly chosen, and afterwards as subsequence within this gene's head domain is randomly selected
    and inverted.

    :param individual: a chromosome
    :return: a tuple of one individual
    """
    if individual.head_length < 2:
        return
    gene = random.choice(individual)
    start, end = _choose_subsequence(gene.head, 2, gene.head_length)
    gene[start: end] = reversed(gene[start: end])
    return individual,


def _choose_donor_donee(individual):
    i1, i2 = random.choices(range(len(individual)), k=2)  # with replacement
    return individual[i1], individual[i2], i1, i2


def isTranspose(individual):
    """
    Perform IS transposition in place

    :param individual: a chromosome
    :return: a tuple of one individual

    An Insertion Sequence (IS) is a random chosen segment across the chromosome, and then the
    IS is copied to be inserted at another position in the head of gene, except the start position. Note that an IS from
    one gene may be copied to a different gene. Typically, `ispb` is set around 0.1.
    """
    donor, donee, i1, i2 = _choose_donor_donee(individual)
    is_start, is_end = _choose_subsequence(donor, 1, donee.head_length - 1)
    is_ = donor[is_start: is_end]
    insertion_pos = random.randint(1, donee.head_length - len(is_))
    donee[:] = donee[:insertion_pos] + is_ + \
        donee[insertion_pos: insertion_pos + donee.head_length - insertion_pos - len(is_)] + donee.tail
    if _DEBUG:
        print('IS transpose: g{}[{}:{}] -> g{}[{}:]'.format(i1, is_start, is_end, i2, insertion_pos))
    return individual,


def risTranspose(individual):
    """
    Perform RIS transposition in place

    :param individual: a chromosome
    :return: a tuple of one individual

    A Root Insertion Sequence (RIS) is a segment of consecutive elements that starts with a
    function. Once an RIS is randomly chosen, it is copied and inserted at the root (first position) of a gene. Note
    that an RIS from one gene may be copied to a different gene. Typically, `rispb` is set around 0.1.
    """
    n_trial = 0
    while n_trial <= 2 * len(individual):
        donor, donee, i1, i2 = _choose_donor_donee(individual)
        # choose a function node randomly to start RIS
        function_indices = [i for i, p in enumerate(donor) if isinstance(p, Function)]
        if not function_indices:  # no functions in this donor, try another
            n_trial += 1
            continue
        ris_start = random.choice(function_indices)
        # determine the length randomly
        length = random.randint(2, min(donee.head_length, len(donor) - ris_start))
        # insert ris at the root of donee
        ris = donor[ris_start: ris_start + length]
        donee[:] = ris + donee[0: donee.head_length - length] + donee.tail
        if _DEBUG:
            print('RIS transpose: g{}[{}:{}] -> g{}[0:]'.format(i1, ris_start, ris_start + length, i2))
        return individual,
    return individual,


def geneTranspose(individual):
    """
    Perform gene transposition in place

    :param individual: a chromosome
    :return: a tuple of one individual

    An entire gene is selected randomly and exchanged with the first gene in the chromosome.
    Obviously, this operation only makes sense if the chromosome is a multigenic one. Typically, `gpb` is set
    around 0.1.

    """
    if len(individual) <= 1:
        return individual,
    source = random.randint(1, len(individual) - 1)
    individual[0], individual[source] = individual[source], individual[0]
    if _DEBUG:
        print('Gene transpose: g0 <-> g{}'.format(source))
    return individual,

