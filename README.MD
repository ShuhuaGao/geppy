# *geppy:* gene expression programming in Python
<img src="docs/source/_static/geppy-icon.png" alt="drawing" width="300px"/>

*geppy* is a computational framework dedicated to [Gene Expression Programming](https://en.wikipedia.org/wiki/Gene_expression_programming) (GEP) in Python,  which is proposed by C. Ferreira  in 2001 [1].  *geppy* is developed in Python 3.

## What is GEP?
Gene Expression Programming (GEP) is a popular and established evolutionary algorithm for automatic generation of computer programs and mathematical models.  It has found wide applications in symbolic regression, classification, automatic model design, combinatorial optimization and real parameter optimization problems [2].

GEP can be seen as a variant of [genetic programming](https://en.wikipedia.org/wiki/Genetic_programming) (GP) and it uses simple linear chromosomes of fixed lengths to encode the genetic information. Though the chromosome (genes) is of fixed length, it can produce expression trees of various sizes thanks to its genotype-phenotype system. Many experiments show that GEP is more efficient than GP, and the trees generated evolved GEP tend to have a smaller size than the ones of GP. 

## *geppy* and [DEAP](https://github.com/DEAP/deap)
*geppy* is built on top of the excellent evolutionary computation framework [DEAP](https://github.com/DEAP/deap) for rapid prototyping and testing of ideas with GEP. DEAP provideds fundamental support for GP, while lacking support for GEP. *geppy* tries the best to follow the style of DEAP and attempts to maintain compatibility with the infrastructure of DEAP. In other words, *geppy* may be considered as a plugin of DEAP to support GEP. If you are familiar with DEAP, then it is easy to grasp *geppy*.

## Features
- Compatilibity with the [DEAP](https://github.com/DEAP/deap) infrastructure and accessibility to DEAP's functionality including:
  - Multi-objective optimisation
  - Parallelization of the evaluations
  - Hall of Fame of the best individuals that lived in the population
  - Checkpoints that take snapshots of a system regularly
- Core data structures in GEP, including gene, chromosome, expression tree, and K-expression.
- Implementation of common mutation, transposition, inversion and crossover operators in GEP.
- Boilerplate algorithms, including the simple GEP algorithm and symbolic regression algorithm (TODO).
- Support numerical constants inference with a third Dc domain in genes. (TODO)
- Visualization of the expression tree.
- Symbolic simplification of a gene, a chromosome, or a K-expression in postprocessing.
- Examples of different applications  in GEP.

## Installation
Currently, *geppy* is still in its alpha phase. If you want to try it, you can install it from sources.
1. First download or clone this repository
```bash
git clone https://github.com/ShuhuaGao/geppy
```
2. Change into the root directory, and iInstall *geppy*
```bash
cd geppy
pip install .
```
(TODO) Later, I will publish *geppy* to pip. 
## Documentation
Check [*geppy* documentation](http://geppy.readthedocs.io/en/latest/) for a comprehensive introduction of its API and typical usages.

TODO: add more examples and tutorials into the documentation.

## Examples

### Symbolic regression
1. [Boolean model identification](./examples/sr/sr_boolean.py)

TODO: add more Jupyter notebook based examples

## Requirements
- Python 3
- [DEAP](https://github.com/DEAP/deap), which should be installed automatically if you haven't got it when installing *geppy*.
- [optional] To visualize the expression tree using the `geppy.export_expression_tree` method, you need the [graphviz](https://pypi.org/project/graphviz/) module.
- [optional] Since GEP/GP doesn't simply the expressions during evolution, its final result may contain many redundancies, like `x + 5 * (2 * x - x - x) - 1`,  which is simply `x - 1`. You may like to simply the final model evolved by GEP with symbolic computation. The corresponding `geppy.simplify` method depends on the [sympy](http://www.sympy.org/en/index.html) package. 

## Reference
The bible of GEP is definitely Ferreira, C.'s monograph: *Ferreira, C. (2006). Gene expression programming: mathematical modeling by an artificial intelligence (Vol. 21). Springer*.

[1] Ferreira, C. (2001). Gene Expression Programming: a New Adaptive Algorithm for Solving Problems. Complex Systems, 13.
[2] Zhong, J., Feng, L., & Ong, Y. S. (2017). Gene expression programming: a survey. IEEE Computational Intelligence Magazine, 12(3), 54-72.