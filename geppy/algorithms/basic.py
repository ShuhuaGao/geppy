# coding=utf-8
"""
.. moduleauthor:: Shuhua Gao

This :mod:`basic` module provides fundamental boilerplate GEP algorithm implementations. After registering proper
operations into a :class:`deap.base.Toolbox` object, the GEP evolution can be simply launched using the present
algorithms. Of course, for complicated problems, you may want to define your own algorithms, and the implementation here
can be used as a reference.
"""
from deap import tools
import random


def _apply_modification(population, operator, pb):
    """
    Apply the modification given by *operator* to each individual in *population* with probability *pb* in place.
    """
    for i in range(len(population)):
        if random.random() < pb:
            population[i], = operator(population[i])
            del population[i].fitness.values


def _apply_crossover(population, operator, pb):
    """
    Mate the *population* in place using *operator* with probability *pb*.
    """
    for i in range(1, len(population), 2):
        if random.random() < pb:
            population[i - 1], population[i] = operator(population[i - 1], population[i])
            del population[i - 1].fitness.values
            del population[i].fitness.values


def _modify(population, toolbox, mutpb, invpb, ispb, rispb, gpb,
            dc_mutpb=0, dc_invpb=0, dc_transpb=0, rnc_mutpb=0, rep=True):
    """
    Modify the population in place and the modified population is also returned.
    """
    # replication, after which all individuals are independent
    if rep:
        offspring = [toolbox.clone(ind) for ind in population]
    else:
        offspring = population

    # mutate
    if hasattr(toolbox, 'mutate'):
        _apply_modification(offspring, toolbox.mutate, mutpb)
    # inversion
    if hasattr(toolbox, 'invert'):
        _apply_modification(offspring, toolbox.invert, invpb)
    # transposition
    if hasattr(toolbox, 'is_transpose'):
        _apply_modification(offspring, toolbox.is_transpose, ispb)
    if hasattr(toolbox, 'ris_transpose'):
        _apply_modification(offspring, toolbox.ris_transpose, rispb)
    if hasattr(toolbox, 'gene_transpose'):
        _apply_modification(offspring, toolbox.gene_transpose, gpb)
    # dc specific
    if hasattr(toolbox, 'mutate_dc'):
        _apply_modification(offspring, toolbox.mutate_dc, dc_mutpb)
    if hasattr(toolbox, 'invert_dc'):
        _apply_modification(offspring, toolbox.invert_dc, dc_invpb)
    if hasattr(toolbox, 'transpose_dc'):
        _apply_modification(offspring, toolbox.transpose_dc, dc_transpb)
    if hasattr(toolbox, 'mutate_rnc_array'):
        _apply_modification(offspring, toolbox.mutate_rnc_array, rnc_mutpb)

    return offspring


def _crossover(population, toolbox, cx1pb, cx2pb, cxgpb, rep=True):
    """
    Perform mating and the children are returned. Every pair of parents will generate two children.
    """
    if rep:
        offspring = [toolbox.clone(ind) for ind in population]
    else:
        offspring = population
    if hasattr(toolbox, 'crossover_one_point'):
        _apply_crossover(offspring, toolbox.crossover_one_point, cx1pb)
    if hasattr(toolbox, 'crossover_two_point'):
        _apply_crossover(offspring, toolbox.crossover_two_point, cx2pb)
    if hasattr(toolbox, 'crossover_gene'):
        _apply_crossover(offspring, toolbox.crossover_gene, cxgpb)
    return offspring


def _local_search(population, toolbox, lspb, rep=True):
    """
    Perform local search in place
    """
    if not hasattr(toolbox, 'local_search'):
        return population
    if rep:
        offspring = [toolbox.clone(ind) for ind in population]
    else:
        offspring = population
    for i in range(len(offspring)):
        if random.random() < lspb:
            offspring[i], = toolbox.local_search(offspring[i])
    return offspring


def gep_simple(population, toolbox, mutation_pb=1, inversion_pb=0.1, is_transposition_pb=0.1, ris_transposition_pb=0.1,
               gene_transposition_pb=0.1, crossover_1p_pb=0.4, crossover_2p_pb=0.2, crossover_gene_pb=0.1,
               n_elites=1, n_generations=100, stats=None, hall_of_fame=None, verbose=__debug__):
    """
    This algorithm performs the simplest gene expression programming. The flowchart of this algorithm can be found
    `here <https://www.gepsoft.com/gxpt4kb/Chapter06/Section1/SS1.htm>`_.
    Refer to Chapter 3 of [FC2006]_ to learn more about this basic algorithm.

    :param population: A list of individuals.
    :param toolbox: A :class:`deap.base.Toolbox` that contains the evolution
                    operators.
    :param mutation_pb: probability of mutating an individual
    :param inversion_pb: probability of inversion
    :param is_transposition_pb: probability of IS transposition
    :param ris_transposition_pb: probability of RIS transposition
    :param gene_transposition_pb: probability of gene transposition
    :param crossover_1p_pb: probability of one-point crossover
    :param crossover_2p_pb: probability of two-point crossover
    :param crossover_gene_pb: probability of gene crossover
    :param n_generations: number of generation.
    :param n_elites: number of elites to be cloned to next generation
    :param stats: a :class:`~deap.tools.Statistics` object that is updated
                  inplace, optional.
    :param hall_of_fame: a :class:`~deap.tools.HallOfFame` object that will
                       contain the best individuals, optional.
    :param verbose: whether or not to log the statistics.
    :returns: The final population
    :returns: A :class:`~deap.tools.Logbook` recording the statistics of the
              evolution process

    .. note::
        This function expects the following aliases to be registered in the toolbox: ``mutate``, ``invert``,
        ``is_transpose``, ``ris_transpose`` and ``gene_transpose`` for modification, ``crossover_one_point``,
        ``crossover_two_point`` and ``crossover_gene`` for crossover, and ``select`` and ``evaluate`` for selection
        and evaluation.

        If an alias is missing in *toolbox*, then it is equivalent to setting a zero probability for that operator.
    """
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    for gen in range(n_generations + 1):

        # evaluate: only evaluate the invalid ones, i.e., no need to reevaluate the unchanged ones
        invalid_individuals = [ind for ind in population if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_individuals)
        for ind, fit in zip(invalid_individuals, fitnesses):
            ind.fitness.values = fit

        # record statistics and log
        if hall_of_fame is not None:
            hall_of_fame.update(population)
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_individuals), **record)
        if verbose:
            print(logbook.stream)

        if gen == n_generations:
            break

        # selection with elitism
        elites = tools.selBest(population, k=n_elites)
        offspring = toolbox.select(population, len(population) - n_elites)

        # modify the offspring after replication
        offspring = _modify(offspring, toolbox, mutation_pb, inversion_pb, is_transposition_pb, ris_transposition_pb,
                            gene_transposition_pb, rep=True)

        # mating
        offspring = _crossover(offspring, toolbox, crossover_1p_pb, crossover_2p_pb, crossover_gene_pb, rep=False)

        # replace the current population with the offsprings
        population = elites + offspring

    return population, logbook


def gep_rnc(population, toolbox, mutation_pb=1, inversion_pb=0.1, is_transposition_pb=0.1, ris_transposition_pb=0.1,
            gene_transposition_pb=0.1, crossover_1p_pb=0.4, crossover_2p_pb=0.2, crossover_gene_pb=0.1,
            dc_mutation_pb=1, dc_inversion_pb=0.1, dc_transposition_pb=0.1, dc_rnc_array_mutation_pb=1,
            n_elites=1, n_generations=100, stats=None, hall_of_fame=None, verbose=__debug__):
    """
    GEP-RNC algorithm, which can evolve random numeric constants (RNCs) effectively. If numerical constants are involved
    in the model, the :class:`~geppy.core.entity.GeneDc` genes should be used accompanied with this algorithm.
    Compared with the basic :func:`gep_simple`
    algorithm, specific operators are used to evolve the Dc domain of :class:`~geppy.core.entity.GeneDc` genes including
    Dc-specific mutation/inversion/transposition and direct mutation of the RNC array associated with each gene. Thus,
    the probability of each of these operators should be specificed in the arguments, which are prefixed with `dc_`.

    Refer to Chapter 5 of [FC2006]_ for more knowledge about GEP-RNC.

    :param population: A list of individuals.
    :param toolbox: A :class:`deap.base.Toolbox` that contains the evolution
                    operators.
    :param mutation_pb: probability of mutating an individual
    :param inversion_pb: probability of inversion
    :param is_transposition_pb: probability of IS transposition
    :param ris_transposition_pb: probability of RIS transposition
    :param gene_transposition_pb: probability of gene transposition
    :param crossover_1p_pb: probability of one-point crossover
    :param crossover_2p_pb: probability of two-point crossover
    :param crossover_gene_pb: probability of gene crossover
    :param dc_mutation_pb: probability of Dc-specific mutation
    :param dc_inversion_pb: probability of Dc-specific inversion
    :param dc_transposition_pb: probability of Dc-specific transposition
    :param dc_rnc_array_mutation_pb: probability of direction mutation of RNC arrays
    :param n_generations: number of generation
    :param n_elites: number of elites to be cloned to next generation
    :param stats: a :class:`~deap.tools.Statistics` object that is updated
                  inplace, optional.
    :param hall_of_fame: a :class:`~deap.tools.HallOfFame` object that will
                       contain the best individuals, optional.
    :param verbose: whether or not to log the statistics.
    :returns: The final population
    :returns: A :class:`~deap.tools.Logbook` recording the statistics of the
              evolution process

    .. note::
        This function expects the following aliases to be registered in the toolbox: ``mutate``, ``invert``,
        ``is_transpose``, ``ris_transpose`` and ``gene_transpose`` for modification, ``crossover_one_point``,
        ``crossover_two_point`` and ``crossover_gene`` for crossover, and ``select`` and ``evaluate`` for selection
        and evaluation. Besides, the four Dc-specific operators are also expected: ``mutate_dc``, ``invert_dc``,
        ``transpose_dc`` and ``mutate_rnc_array_dc``. Besides, ``toolbox.local_search`` is used to perform local search
        for RNC.

        If an alias is missing in *toolbox*, then it is equivalent to setting a zero probability for that operator.
    """
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    for gen in range(n_generations + 1):

        # evaluate: only evaluate the invalid ones, i.e., no need to reevaluate the unchanged ones
        invalid_individuals = [ind for ind in population if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_individuals)
        for ind, fit in zip(invalid_individuals, fitnesses):
            ind.fitness.values = fit

        # record statistics and log
        if hall_of_fame is not None:
            hall_of_fame.update(population)
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_individuals), **record)
        if verbose:
            print(logbook.stream)

        if gen == n_generations:
            break

        # selection with elitism
        elites = tools.selBest(population, k=n_elites)
        offspring = toolbox.select(population, len(population) - n_elites)

        # modify the offspring after replication
        offspring = _modify(offspring, toolbox, mutation_pb, inversion_pb, is_transposition_pb, ris_transposition_pb,
                            gene_transposition_pb, dc_mutation_pb, dc_inversion_pb, dc_transposition_pb,
                            dc_rnc_array_mutation_pb, rep=True)

        # mating
        offspring = _crossover(offspring, toolbox, crossover_1p_pb, crossover_2p_pb, crossover_gene_pb, rep=False)

        # replace the current population with the offsprings
        population = elites + offspring

    return population, logbook
