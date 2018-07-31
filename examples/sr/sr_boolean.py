# coding=utf-8
"""
Symbolic regression for Boolean expressions.
"""
import geppy as gep
from deap import base, creator, tools
import operator
import numpy
import itertools
import random


# the true model
def f(a, b, c, d):
    return (a and d) or not c


# generate the training set which contains all the 16 samples
X = []
Y = []
for a, b, c, d in itertools.product([True, False], repeat=4):
    X.append((a, b, c, d))
    Y.append(f(a, b, c, d))


# define the primitive set
pset = gep.PrimitiveSet('Main', input_names=['a', 'b', 'c', 'd'])
pset.add_function(operator.and_, 2)
pset.add_function(operator.or_, 2)
pset.add_function(operator.not_, 1)


# Individual class, which is a subclass of gep.Chromosome
# We want to minimize the fitting error. Check http://deap.readthedocs.io/en/master/api/base.html#deap.base.Fitness.
creator.create("FitnessMax", base.Fitness, weights=(1,))
creator.create("Individual", gep.Chromosome, fitness=creator.FitnessMax)

# toolbox: register necessary operations to create individuals and population
h = 7   # head length
n_genes = 2
toolbox = base.Toolbox()
toolbox.register('gene_gen', gep.Gene, pset=pset, head_length=h)
toolbox.register('individual', creator.Individual, gene_gen=toolbox.gene_gen, n_genes=n_genes, linker=operator.or_)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register('compile', gep.compile_, pset=pset)


# fitness evaluation function
def evaluate(individual):
    """
    Compute the sum of absolute errors (SAE).
    """
    func = toolbox.compile(individual)  # the function evolved by GEP
    n_correct = 0
    for (a, b, c, d), y in zip(X, Y):
        prediction = func(a, b, c, d)
        if prediction == y:
            n_correct += 1
    return n_correct,


toolbox.register('evaluate', evaluate)
# selection, modification and crossover operators
toolbox.register('select', tools.selRoulette)
toolbox.register('mutate', gep.mutate_uniform, pset=pset, indpb=2 / (2 * h + 1))
toolbox.register('invert', gep.invert)
toolbox.register('isTranspose', gep.is_transpose)
toolbox.register('risTranspose', gep.ris_transpose)
toolbox.register('geneTranspose', gep.gene_transpose)
toolbox.register('mate1p', gep.crossover_one_point)
toolbox.register('mate2p', gep.crossover_two_point)
toolbox.register('mateg', gep.crossover_gene)


def main():
    # random.seed(2)  # for reproduction purpose
    n_pop = 50
    n_gen = 50

    pop = toolbox.population(n=n_pop)
    hof = tools.HallOfFame(3)

    stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)

    # start evolution
    pop, log = gep.gep_simple(pop, toolbox, mutpb=0.9, invpb=0.1, ispb=0.1, rispb=0.1, gpb=0.1,
                              cx1pb=0.4, cx2pb=0.2, cxgpb=0.1,
                              n_gen=n_gen, n_elites=2,
                              stats=stats, halloffame=hof)
    return pop, log, hof


if __name__ == '__main__':
    pop, log, hof = main()
    rename_labels = {'and_': '&', 'or_': '|', 'not_': '~'}
    nodes, edges, labels = gep.graph(hof[0], rename_labels)
    gep.export_expression_tree(hof[0], rename_labels, file='./data/sr_boolean_tree.png')

    best_individual = hof[0]
    solution = gep.simplify(hof[0])
    print(solution)