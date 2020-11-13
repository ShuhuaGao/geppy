# coding=utf-8
"""
.. moduleauthor:: Shuhua Gao

The module :mod:`mutation` provides mutation related genetic modification operators in GEP, including point-mutation,
transposition and inversion. Please refer to Chapter 3 of [FC2006]_ for more details.

.. note::
    1. Operators that can be applied to both :class:`~geppy.core.entity.Gene` and :class:`~geppy.core.entity.GeneDc`:
       :func:`mutate_uniform`, :func:`invert`, :func:`is_transpose`, :func:`ris_transpose`, :func:`gene_transpose`
    2. Operators that end with '_dc' are specially designed to handle the Dc domain in
       :class:`~geppy.core.entity.GeneDc`.
"""
import random
from ..core.symbol import Function, EphemeralTerminal
from ._util import _choose_a_terminal

_DEBUG = False


def _choose_function(pset):
    return random.choice(pset.functions)


def _choose_terminal(pset):
    return _choose_a_terminal(pset.terminals)


def _choose_subsequence(seq, min_length=1, max_length=-1):
    if max_length <= 0:
        max_length = len(seq)
    length = random.randint(min_length, max_length)
    start = random.randint(0, len(seq) - length)
    return start, start + length


def _choose_subsequence_indices(i, j, min_length=1, max_length=-1):
    """
    Choose a subsequence from [i, j] (both included) and return the subsequence boundaries [a, b] (both included).
    Additionally, the min_length and max length of [a, b], i.e., b - a + 1 can be specified.
    """
    if max_length <= 0:
        max_length = j - i + 1
    length = random.randint(min_length, max_length)
    start = random.randint(i, j - length + 1)
    return start, start + length - 1


def mutate_uniform(individual, pset, ind_pb='2p'):
    """
    Uniform point mutation. For each symbol (primitive) in *individual*, change it to another randomly chosen symbol
    from *pset* with the probability *indpb*. A symbol may be a function or a terminal.

    :param individual: :class:`~geppy.core.entity.Chromosome`, the chromosome to be mutated.
    :param pset: :class:`~geppy.core.entity.PrimitiveSet`, a primitive set
    :param ind_pb: float or str, default '2p'. Probability of mutating each symbol.
        If *ind_pb* is given as a string ending with 'p',
        then it indicates the expected number of point mutations among all the symbols in *individual*. For example,
        if the total length of each gene of *individual* is `l` and there are `m` genes in total, then by passing
        `ind_pb='1.5p'` we specify approximately `ind_pb=1.5/(l*m)`.
    :return: A tuple of one chromosome

    It is typical to set a mutation rate *indpb* equivalent to two one-point mutations per chromosome. That is,
    ``indpb = 2 / len(chromosome) * len(gene)``.
    """
    if isinstance(ind_pb, str):
        assert ind_pb.endswith('p'), "ind_pb must end with 'p' if given in a string form"
        length = individual[0].head_length + individual[0].tail_length
        ind_pb = float(ind_pb.rstrip('p')) / (len(individual) * length)
    for gene in individual:
        # mutate the gene with the associated pset
        # head: any symbol can be changed into a function or a terminal
        for i in range(gene.head_length):
            if random.random() < ind_pb:
                if random.random() < 0.5:  # to a function
                    gene[i] = _choose_function(pset)
                else:
                    gene[i] = _choose_terminal(pset)
        # tail: only change to another terminal
        for i in range(gene.head_length, gene.head_length + gene.tail_length):
            if random.random() < ind_pb:
                gene[i] = _choose_terminal(pset)
    return individual,


def mutate_uniform_dc(individual, ind_pb='1p'):
    """
    Dc-specific mutation. This operator changes one index stored in the Dc domain to another index in place. The indices
    in the Dc domain can later be used to retrieve numerical constants from the gene associated
    :meth:`~geppy.core.entity.GeneDc.rnc_array`.

    :param individual: :class:`~geppy.core.entity.Chromosome`, a chromosome, which contains genes of type
        :class:`~geppy.core.entity.GeneDc`
    :param ind_pb: float or str, default '1p'. Probability of mutating each index/position in the Dc domain. If given
        in str like 'xp', then `x` point mutations for expected for the Dc domains in total of *individual*. For
        instance, if there are `d` Dc elements in the whole chromosome, then `1.5p` is equal to `ind_pb = 1.5 / d`.
    :return: a tuple of one chromosome

    It is typical to set a mutation rate *indpb* equivalent to two one-point mutations per chromosome. That is,
    ``indpb = 2 / len(chromosome) * len(gene)``.
    """
    if isinstance(ind_pb, str):
        assert ind_pb.endswith('p'), "ind_pb must end with 'p' if given in a string form"
        ind_pb = float(ind_pb.rstrip('p')) / (individual[0].dc_length * len(individual))
    # change the position in the Dc domain to another one
    for g in individual:  # g: GeneDc
        start = g.head_length + g.tail_length
        end = start + g.dc_length
        for i in range(start, end):
            if random.random() < ind_pb:
                g[i] = random.randint(0, len(g.rnc_array) - 1)
    return individual,


def invert(individual):
    """
    A gene is randomly chosen, and afterwards a subsequence within this gene's head domain is randomly selected
    and inverted.

    :param individual: :class:`~geppy.core.entity.Chromosome`, a chromosome
    :return: a tuple of one individual

    Typically, a small inversion rate of 0.1 is used.
    """
    if individual.head_length < 2:
        return individual,
    gene = random.choice(individual)
    start, end = _choose_subsequence(gene.head, 2, gene.head_length)
    gene[start: end] = reversed(gene[start: end])
    if _DEBUG:
        print('invert [{}: {}]'.format(start, end))
    return individual,


def invert_dc(individual):
    """
    Dc-specific inversion. A gene is randomly chosen, and afterwards a subsequence within this gene's Dc domain is
    randomly selected and inverted.

    :param individual: :class:`~geppy.core.entity.Chromosome`, a chromosome, which contains genes of type
        :class:`~geppy.core.entity.GeneDc`
    :return: a tuple of one individual

    Typically, a small inversion rate of 0.1 is used.
    """
    g = random.choice(individual)
    if g.dc_length < 2:
        return individual,
    start = g.head_length + g.tail_length
    end = start + g.dc_length - 1
    inv_start, inv_end = _choose_subsequence_indices(start, end, min_length=2)
    g[inv_start: inv_end] = reversed(g[inv_start: inv_end])
    return individual,


def _choose_donor_donee(individual):
    i1, i2 = random.choices(range(len(individual)), k=2)  # with replacement
    return individual[i1], individual[i2], i1, i2


def is_transpose(individual):
    """
    Perform IS transposition in place

    :param individual: :class:`~geppy.core.entity.Chromosome`, a chromosome
    :return: a tuple of one individual

    An Insertion Sequence (IS) is a random chosen segment across the chromosome, and then the
    IS is copied to be inserted at another position in the head of gene, except the start position. Note that an IS from
    one gene may be copied to a different gene.
    Typically, the IS transposition rate is set around 0.1.
    """
    donor, donee, i1, i2 = _choose_donor_donee(individual)
    a, b = _choose_subsequence_indices(0, donor.head_length + donor.tail_length - 1, max_length=donee.head_length - 1)
    is_start, is_end = a, b + 1
    is_ = donor[is_start: is_end]
    insertion_pos = random.randint(1, donee.head_length - len(is_))
    donee[:] = donee[:insertion_pos] + is_ + \
               donee[insertion_pos: insertion_pos + donee.head_length - insertion_pos - len(is_)] + \
               donee[donee.head_length:]
    if _DEBUG:
        print('IS transpose: g{}[{}:{}] -> g{}[{}:]'.format(i1, is_start, is_end, i2, insertion_pos))
    return individual,


def ris_transpose(individual):
    """
    Perform RIS transposition in place

    :param individual: :class:`~geppy.core.entity.Chromosome`, a chromosome
    :return: a tuple of one individual

    A Root Insertion Sequence (RIS) is a segment of consecutive elements that starts with a
    function. Once an RIS is randomly chosen, it is copied and inserted at the root (first position) of a gene. Note
    that an RIS from one gene may be copied to a different gene.
    Typically, the RIS transposition rate is set around 0.1.
    """
    n_trial = 0
    while n_trial <= 2 * len(individual):
        donor, donee, i1, i2 = _choose_donor_donee(individual)
        # choose a function node randomly to start RIS
        function_indices = [i for i, p in enumerate(donor.head) if isinstance(p, Function)]
        if not function_indices:  # no functions in this donor, try another
            n_trial += 1
            continue
        ris_start = random.choice(function_indices)
        # determine the length randomly
        length = random.randint(2, min(donee.head_length, donor.head_length + donor.tail_length - ris_start))
        # insert ris at the root of donee
        ris = donor[ris_start: ris_start + length]
        donee[:] = ris + donee[0: donee.head_length - length] + donee[donee.head_length:]
        if _DEBUG:
            print('RIS transpose: g{}[{}:{}] -> g{}[0:]'.format(i1, ris_start, ris_start + length, i2))
        return individual,
    return individual,


def gene_transpose(individual):
    """
    Perform gene transposition in place

    :param individual: :class:`~geppy.core.entity.Chromosome`, a chromosome
    :return: a tuple of one individual

    An entire gene is selected randomly and exchanged with the first gene in the chromosome.
    Obviously, this operation only makes sense if the chromosome is a multigenic one. Typically, the gene transposition
    rate is set around 0.1.

    """
    if len(individual) <= 1:
        return individual,
    source = random.randint(1, len(individual) - 1)
    individual[0], individual[source] = individual[source], individual[0]
    if _DEBUG:
        print('Gene transpose: g0 <-> g{}'.format(source))
    return individual,


def transpose_dc(individual):
    """
    This operator randomly chooses the chromosome, the gene with its respective Dc to be subjected to transposition,
    the start and termination points of the transposon, and the target site; then it inserts a copy of the transposon
    from the place of origin to the target site, while the old content at the target site are shifted to the right
    and the length of the Dc domain is maintained.

    :param individual: :class:`~geppy.core.entity.Chromosome`, a chromosome, which contains genes of type
        :class:`~geppy.core.entity.GeneDc`
    :return: a tuple of one individual

    Typically, a small transposition rate of 0.1 is used. Note that a transposon from one gene's Dc may be copied and
    inserted into another gene's Dc, which means the genes in a multigenic chromosome *individual*
    must have an identical structure.
    """
    donor, donee, i, j = _choose_donor_donee(individual)
    dc_start = donor.head_length + donor.tail_length
    dc_end = donor.head_length + donor.tail_length + donor.dc_length - 1  # included
    # choose a transposon from the donor Dc
    t_start, t_end = _choose_subsequence_indices(dc_start, dc_end)
    # choose the insertion point in donee Dc
    insertion_pos = random.randint(dc_start, dc_end - (t_end - t_start))
    donee[:] = donee[:insertion_pos] + donor[t_start: t_end + 1] \
               + donee[insertion_pos: insertion_pos + len(donee) - (insertion_pos + t_end - t_start + 1)]
    if _DEBUG:
        print('Dc transpose: g{}[{}:{}] -> g{}[{}] (both included)'.format(i, t_start, t_end, j, insertion_pos))
    return individual,


def mutate_rnc_array_dc(individual, rnc_gen, ind_pb='1p'):
    """
    Direct mutation of RNCs, which changes the values in a gene's RNC array
    :meth:`~geppy.core.entity.GeneDc.rnc_array` randomly.

    :param individual: :class:`~geppy.core.entity.Chromosome`, a chromosome, which contains genes of type
        :class:`~geppy.core.entity.GeneDc`.
    :param rnc_gen: callable, which returns a random numerical constant by calling ``rnc_gen()``.
    :param ind_pb: float or str, default '1p'. Use a float number to specify the probability of each RNC in the array
        being mutated. Alternatively, if a str ending with 'p' is given in the form like 'xp',  where
        'x' is a numerical number, then it specifies how many point mutations are expected for the RNC arrays
        in *individual*. For example, if there are *d* RNCs inside the arrays in total in *individual*,
        then '1.5p' is approximately equal to `ind_pb = 1.5 / d`.
    :return: a tuple of one individual

    The genetic operators :func:`mutate_uniform_dc`, :func:`transpose_dc` and :func:`invert_dc` actually only move the
    random numerical constants around without generating new numerical values. This method :func:`mutate_rnc_array_dc`
    can replace the value of a particular numerical constant by another in the
    :meth:`~geppy.core.entity.GeneDc.rnc_array`.

    Refer to section 5.4.4 of [FC2006]_ for more details.
    """
    if isinstance(ind_pb, str):
        assert ind_pb.endswith('p'), "ind_pb must end with 'p' if given in a string form"
        n_rnc_total = len(individual[0].rnc_array) * len(individual)
        if n_rnc_total == 0:
            return individual,
        ind_pb = float(ind_pb.rstrip('p')) / n_rnc_total

    for g in individual:
        for i in range(len(g.rnc_array)):
            if random.random() < ind_pb:
                g.rnc_array[i] = rnc_gen()
    return individual,


def mutate_uniform_ephemeral(individual, ind_pb='1p'):
    """
    This operator is specially designed for ephemeral numerical constant (ENC) mutation.
    It attempts to change the value of the ephemeral constants in the *individual* if any by calling `e.update_value()`,
    where `e` is an ephemeral inside *individual*. See :meth:`~geppy.core.symbol.EphemeralTerminal.update_value` for details.

    :param individual: :class:`~geppy.core.entity.Chromosome`, a chromosome
    :param ind_pb: float or str, default '1p'. If a float is given, then *ind_pb* specifies the probability of each
        ephemeral in the *individual* being mutated. If a str ending with 'p' is given in the form like 'xp',  where
        'x' is a numerical number, then it
        specifies how many point mutations are expected for the ephemeral constants. For example, if there are *d*
        ephemeral constants present in total in *individual*, then '1.5p' is approximately equal to `ind_pb = 1.5 / d`.

    :return: a tuple of one individual
    """
    ephemerals = [p for g in individual for p in g if isinstance(p, EphemeralTerminal)]
    if len(ephemerals) == 0:
        return individual,
    if isinstance(ind_pb, str):
        assert ind_pb.endswith('p'), "ind_pb must end with 'p' if given in a string form"
        np = float(ind_pb.rstrip('p'))
        ind_pb = np / len(ephemerals)
    for e in ephemerals:
        if random.random() < ind_pb:
            e.update_value()
    return individual,
