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


def _modify(population, toolbox, mutpb, invpb, ispb, rispb, gpb, rep=True):
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
    if hasattr(toolbox, 'isTranspose'):
        _apply_modification(offspring, toolbox.isTranspose, ispb)
    if hasattr(toolbox, 'risTranspose'):
        _apply_modification(offspring, toolbox.risTranspose, rispb)
    if hasattr(toolbox, 'geneTranspose'):
        _apply_modification(offspring, toolbox.geneTranspose, gpb)

    return offspring


def _crossover(population, toolbox, cx1pb, cx2pb, cxgpb, rep=True):
    """
    Perform mating and the children are returned. Every pair of parents will generate two children.
    """
    if rep:
        offspring = [toolbox.clone(ind) for ind in population]
    else:
        offspring = population
    _apply_crossover(offspring, toolbox.mate1p, cx1pb)
    _apply_crossover(offspring, toolbox.mate2p, cx2pb)
    _apply_crossover(offspring, toolbox.mateg, cxgpb)
    return offspring


def gep_simple(population, toolbox, mutpb, invpb, ispb, rispb, gpb, cx1pb, cx2pb, cxgpb,
               n_elites, n_gen, stats=None, halloffame=None, verbose=__debug__):
    """
    This sr and the simplest gene expression algorithm. The flowchart of this algorithm can be found
    `here <https://www.gepsoft.com/gxpt4kb/Chapter06/Section1/SS1.htm>`_.

    :param population: A list of individuals.
    :param toolbox: A :class:`deap.base.Toolbox` that contains the evolution
                    operators.
    :param mutpb: probability of mutating an individual
    :param invpb: probability of inversion
    :param ispb: probability of IS transposition
    :param rispb: probability of RIS transposition
    :param gpb: probability of gene transposition
    :param cx1pb: probability of one-point crossover
    :param cx2pb: probability of two-point crossover
    :param cxgpb: probability of gene crossover
    :param n_gen: number of generation.
    :param n_elites: number of elites to be cloned to next generation
    :param stats: a :class:`~deap.tools.Statistics` object that is updated
                  inplace, optional.
    :param halloffame: a :class:`~deap.tools.HallOfFame` object that will
                       contain the best individuals, optional.
    :param verbose: whether or not to log the statistics.
    :returns: The final population
    :returns: A :class:`~deap.tools.Logbook` recording the statistics of the
              evolution process

    .. note::
        This function expects the following aliases to be registered in the toolbox: :meth:`toolbox.mutate`,
        :meth:`toolbox.invert`, :meth:`toolbox.isTranspose`, :meth:`toolbox.risTranspose`,
        :meth:`toolbox.geneTranspose` :meth:`toolbox.mate1p`, :meth:`toolbox.mate2p`, :meth:`toolbox.mateg`
        for mutation and crossover, and
        :meth:`toolbox.select`, :meth:`toolbox.evaluate` for selection and evaluation. If an alias is missing
        in *toolbox*, then it is equivalent to setting a zero probability.
    """
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    for gen in range(n_gen + 1):

        # evaluate: only evaluate the invalid ones, i.e., no need to reevaluate the unchanged ones
        invalid_individuals = [ind for ind in population if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_individuals)
        for ind, fit in zip(invalid_individuals, fitnesses):
            ind.fitness.values = fit

        # record statistics and log
        if halloffame is not None:
            halloffame.update(population)
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_individuals), **record)
        if verbose:
            print(logbook.stream)

        if gen == n_gen:
            break

        # selection with elitism
        elites = tools.selBest(population, k=n_elites)
        offspring = toolbox.select(population, len(population) - n_elites)

        # modify the offspring after replication
        offspring = _modify(offspring, toolbox, mutpb, invpb, ispb, rispb, gpb, rep=True)

        # mating
        offspring = _crossover(offspring, toolbox, cx1pb, cx2pb, cxgpb, rep=False)

        # replace the current population with the offsprings
        population = elites + offspring

    return population, logbook




