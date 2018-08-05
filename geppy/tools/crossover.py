# coding=utf-8
"""
.. moduleauthor:: Shuhua Gao

The module :mod:`crossover` provides crossover (mating, recombination) related genetic modifications in GEP,
including one-point and two-point crossover, and gene crossover between multigenic chromosomes.
Please refer to Chapter 3 of [FC2006]_ for more details.

.. note::
    All the recombination operators including :func:`crossover_one_point`, :func:`crossover_two_point`, and
    :func:`crossover_gene` can be applied to both :class:`~geppy.core.entity.Gene`
    and :class:`~geppy.core.entity.GeneDc`.

.. [FC2006] Ferreira, Cândida. Gene expression programming: mathematical modeling by an artificial
    intelligence. Vol. 21. Springer, 2006.

"""
import random

_DEBUG = False


def crossover_one_point(ind1, ind2):
    """
    Execute one-point recombination of two individuals. The two individuals are modified in place, and the two children
    are returned.

    :param ind1: The first individual (chromosome) participating in the crossover.
    :param ind2: The second individual (chromosome) participating in the crossover.
    :return: A tuple of two children individuals.

    Note the crossover can happen at any point across the whole chromosome and thus entire genes may be also exchanged
    between the two parents if they are multigenic chromosomes.
    """
    assert len(ind1) == len(ind2)
    # the gene containing the recombination point, and the point index in the gene
    which_gene = random.randint(0, len(ind1) - 1)
    which_point = random.randint(0, len(ind1[which_gene]) - 1)
    # exchange the upstream materials
    ind1[:which_gene], ind2[:which_gene] = ind2[:which_gene], ind1[:which_gene]
    ind1[which_gene][:which_point + 1], ind2[which_gene][:which_point + 1] = \
        ind2[which_gene][:which_point + 1], ind1[which_gene][:which_point + 1]
    if _DEBUG:
        print('cxOnePoint: g{}[{}]'.format(which_gene, which_point))
    return ind1, ind2


def crossover_two_point(ind1, ind2):
    """
    Execute two-point recombination of two individuals. The two individuals are modified in place, and the two children
    are returned. The materials between two randomly chosen points are swapped to generate two children.

    :param ind1: The first individual (chromosome) participating in the crossover.
    :param ind2: The second individual (chromosome) participating in the crossover.
    :return: A tuple of two individuals.

    Note the crossover can happen at any point across the whole chromosome and thus entire genes may be also exchanged
    between the two parents if they are multigenic chromosomes.
    """
    assert len(ind1) == len(ind2)
    # the two genes containing the two recombination points
    g1, g2 = random.choices(range(len(ind1)), k=2)  # with replacement, thus g1 may be equal to g2
    if g2 < g1:
        g1, g2 = g2, g1
    # the two points in g1 and g2
    p1 = random.randint(0, len(ind1[g1]) - 1)
    p2 = random.randint(0, len(ind1[g2]) - 1)
    # change the materials between g1->p1 and g2->p2: first exchange entire genes, then change partial genes at g1, g2
    if g1 == g2:
        if p1 > p2:
            p1, p2 = p2, p1
        ind1[g1][p1: p2+1], ind2[g2][p1: p2+1] = ind2[g2][p1: p2+1], ind1[g1][p1: p2+1]
    else:
        ind1[g1 + 1: g2], ind2[g1 + 1: g2] = ind2[g1 + 1: g2], ind1[g1 + 1: g2]
        ind1[g1][p1:], ind2[g1][p1:] = ind2[g1][p1:], ind1[g1][p1:]
        ind1[g2][:p2 + 1], ind2[g2][:p2 + 1] = ind2[g2][:p2 + 1], ind1[g2][:p2 + 1]
    if _DEBUG:
        print('cxTwoPoint: g{}[{}], g{}[{}]'.format(g1, p1, g2, p2))
    return ind1, ind2


def crossover_gene(ind1, ind2):
    """
    Entire genes are exchanged between two parent chromosomes. The two individuals are modified in place, and the two
    children are returned.

    :param ind1: The first individual (chromosome) participating in the crossover.
    :param ind2: The second individual (chromosome) participating in the crossover.
    :return: a tuple of two children individuals

    This operation has no effect if the chromosome has only one gene. Typically, a gene recombination rate
    around 0.2 is used.
    """
    assert len(ind1) == len(ind2)
    pos1, pos2 = random.choices(range(len(ind1)), k=2)
    ind1[pos1], ind2[pos2] = ind2[pos2], ind1[pos1]
    if _DEBUG:
        print('cxGene: ind1[{}] <--> ind2[{}]'.format(pos1, pos2))
    return ind1, ind2
