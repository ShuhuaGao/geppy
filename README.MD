# *geppy:* a gene expression programming framework in Python
<img src="https://github.com/ShuhuaGao/geppy/blob/master/docs/source/_static/geppy-icon.png" width="300">

*geppy* is a computational framework dedicated to [Gene Expression Programming](https://en.wikipedia.org/wiki/Gene_expression_programming) (GEP),  which is proposed by C. Ferreira  in 2001 [1].  *geppy* is developed in Python 3.

## What is GEP?
Gene Expression Programming (GEP) is a popular and established evolutionary algorithm for automatic generation of computer programs and mathematical models.  It has found wide applications in symbolic regression, classification, automatic model design, combinatorial optimization and real parameter optimization problems [2].

GEP can be seen as a variant of the traditional  [genetic programming](https://en.wikipedia.org/wiki/Genetic_programming) (GP) and it uses simple linear chromosomes of fixed lengths to encode the genetic information. Though the chromosome (genes) is of fixed length, it can produce expression trees of various sizes thanks to its genotype-phenotype expressio system. Many experiments show that GEP is more efficient than GP, and the trees evolved by GEP tend to have a smaller size than the ones of GP. 

## *geppy* and [DEAP](https://github.com/DEAP/deap)
*geppy* is built on top of the excellent evolutionary computation framework [DEAP](https://github.com/DEAP/deap) for rapid prototyping and testing of ideas with GEP. DEAP provides fundamental support for GP, while lacking support for GEP. *geppy* tries the best to follow the style of DEAP and attempts to maintain compatibility with the major infrastructure of DEAP. In other words, to some degree *geppy* may be considered as a plugin of DEAP to specially support GEP. If you are familiar with DEAP, then it is easy to grasp *geppy*. Besides, comprehensive [documentation](https://geppy.readthedocs.io/en/latest/) is also available.

## Features
- Compatibility with the [DEAP](https://github.com/DEAP/deap) infrastructure and easy accessibility to DEAP's functionality including:
  - Multi-objective optimisation
  - Straightforward parallelization of fitness evaluations for speedup
  - Hall of Fame of the best individuals that lived in the population
  - Checkpoints that take snapshots of a system regularly
  - Statistics and logging
- Core data structures in GEP, including the gene, chromosome, expression tree, and K-expression.
- Implementation of common mutation, transposition, inversion and crossover operators in GEP.
- Boilerplate algorithms, including  the standard GEP algorithm and advanced algorithms integrating a local optimizer for numerical constant optimization.
- Support numerical constants inference with a third Dc domain in genes: the GEP-RNC algorithm.
- Flexible built-in algorithm interface, which can support an arbitrary number of custom mutation and crossover-like operators.
- Visualization of the expression tree.
- Symbolic simplification of a gene, a chromosome, or a K-expression in postprocessing.
- Examples of different applications using  GEP with detailed comments in Jupyter notebook.

## Installation
### From PyPI (recommended)
```bash
pip install geppy
```
### From source
You can install it from sources.
1. First download or clone this repository
```bash
git clone https://github.com/ShuhuaGao/geppy
```
2. Change into the root directory, i.e., the one containing the *setup.py* file and install *geppy* using *pip*
```bash
cd geppy
pip install .
```
## Documentation
Check [*geppy* documentation](https://geppy.readthedocs.io/en/latest/) for GEP theory and tutorials as well as a comprehensive introduction of *geppy*'s API and typical usages with comprehensive tutorials and examples.

## Examples
A getting started example is presented in the Jupyter notebook [Boolean model identification](./examples/sr/Boolean_function_identification.ipynb), which infers a Boolean function from given input-output data with GEP. More examples are listed in the following.
### Simple symbolic regression
1. [Boolean model identification](./examples/sr/Boolean_function_identification.ipynb) (Getting started with no constants involved)
2. [Simple mathematical expression inference](./examples/sr/numerical_expression_inference-ENC.ipynb) (Constants finding with ephemeral random constants (ERC))
3. [Simple mathematical expression inference with the GEP-RNC algorithm](./examples/sr/numerical_expression_inference-RNC.ipynb) (Demonstrating the GEP-RNC algorithm for numerical constant evolution)
### Advanced symbolic regression
1. [Improving symbolic regression with linear scaling](./examples/sr/numerical_expression_inference-Linear_scaling.ipynb) (Use the linear scaling technique to evolve models with continuous real constants more efficiently)

2. [Use the GEP-RNC algorithm with linear scaling on  the UCI Power Plant dataset](https://github.com/ShuhuaGao/geppy/blob/master/examples/sr/GEP_RNC_for_ML_with_UCI_Power_Plant_dataset.ipynb) See how to apply GEP based symbolic regression on a real machine learning dataset.


## Requirements
- Python 3.6 and afterwards
- [DEAP](https://github.com/DEAP/deap), which should be installed automatically if you haven't got it when installing *geppy*.
- [optional] To visualize the expression tree using the `geppy.export_expression_tree` method, you need the [graphviz](https://pypi.org/project/graphviz/) module.
- [optional] Since GEP/GP doesn't simplify the expressions during evolution, its final result may contain many redundancies, and the tree can be very large,  like `x + 5 * (2 * x - x - x) - 1`,  which is simply `x - 1`. You may like to simplify the final model evolved by GEP with symbolic computation to get better understanding of this model. The corresponding `geppy.simplify` method depends on the [sympy](http://www.sympy.org/en/index.html) package. 

## Common pitfalls in using GP
Always keep in mind that evolution is random. Thus, any values may be input into a function. If issues like "overflow", "nan", or "not a number", or unreasonally huge values are encounterred, the most possible reason is that you did not protect a possibly _dangerous_ function. For example, if the `sqrt` function lies in the function set, then in evaluating one individual evolved by `geppy` (or GP in general), it is likely that a negative input `sqrt(-1.24)` may happen. 

Refer to issues [#28](https://github.com/ShuhuaGao/geppy/issues/28) [#26](https://github.com/ShuhuaGao/geppy/issues/26) [#4](https://github.com/ShuhuaGao/geppy/issues/4) for more details.

## Reference
The bible of GEP is definitely Ferreira, C.'s monograph: **Ferreira, C. (2006). Gene expression programming: mathematical modeling by an artificial intelligence (Vol. 21). Springer**.

You can also get a lot of papers/documents by Googling 'gene expression programming'.

[1] Ferreira, C. (2001). Gene Expression Programming: a New Adaptive Algorithm for Solving Problems. Complex Systems, 13.
[2] Zhong, J., Feng, L., & Ong, Y. S. (2017). Gene expression programming: a survey. IEEE Computational Intelligence Magazine, 12(3), 54-72.

## How to cite *geppy*
If you find *geppy* useful in your projects, please cite it such that more researchers/engineers will know it. A BibTeX entry for *geppy* is given below.
```
@misc{geppy_2020,
    author       = {Shuhua Gao},
    title        = {{geppy: a Python framework for gene expression programming }},
    month        = July,
    year         = 2020,
    doi          = {10.5281/zenodo.3946297},
    version      = {0.1},
    publisher    = {Zenodo},
    url          = {https://github.com/ShuhuaGao/geppy}
    }
```
Alternatively, if you want a more *academic* citation, you may cite our relevant paper
```
@ARTICLE{learn_async,
  author={S. {Gao} and C. {Sun} and C. {Xiang} and K. {Qin} and T. H. {Lee}},
  journal={IEEE Transactions on Cybernetics}, 
  title={Learning Asynchronous Boolean Networks From Single-Cell Data Using Multiobjective Cooperative Genetic Programming}, 
  year={2020},
  volume={},
  number={},
  pages={1-15},
  doi={10.1109/TCYB.2020.3022430}}
```
