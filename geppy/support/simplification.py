# coding=utf-8
"""
.. moduleauthor:: Shuhua Gao

This module :mod:`simplification` provides utility functions for symbolic simplification of GEP individuals, which may
be used in postprocessing.
"""
import math
import operator
from ..core.entity import KExpression, Chromosome, Gene
from ..core.symbol import Function, Terminal, SymbolTerminal

import sympy as sp


DEFAULT_SYMBOLIC_FUNCTION_MAP = {
    operator.and_.__name__: sp.And,
    operator.or_.__name__: sp.Or,
    operator.not_.__name__: sp.Not,
    operator.add.__name__: operator.add,
    operator.sub.__name__: operator.sub,
    operator.mul.__name__: operator.mul,
    operator.neg.__name__: operator.neg,
    operator.pow.__name__: operator.pow,
    operator.abs.__name__: operator.abs,
    operator.floordiv.__name__: operator.floordiv,
    operator.truediv.__name__: operator.truediv,
    'protected_div': operator.truediv,
    math.log.__name__: sp.log,
    math.sin.__name__: sp.sin,
    math.cos.__name__: sp.cos,
    math.tan.__name__: sp.tan
}
"""
Currently, it is defined as::

    DEFAULT_SYMBOLIC_FUNCTION_MAP = {
        operator.and_.__name__: sp.And,
        operator.or_.__name__: sp.Or,
        operator.not_.__name__: sp.Not,
        operator.add.__name__: operator.add,
        operator.sub.__name__: operator.sub,
        operator.mul.__name__: operator.mul,
        operator.neg.__name__: operator.neg,
        operator.pow.__name__: operator.pow,
        operator.abs.__name__: operator.abs,
        operator.floordiv.__name__: operator.floordiv,
        operator.truediv.__name__: operator.truediv,
        'protected_div': operator.truediv,
        math.log.__name__: sp.log,
        math.sin.__name__: sp.sin,
        math.cos.__name__: sp.cos,
        math.tan.__name__: sp.tan
    }
"""


def _simplify_kexpression(expr, symbolic_function_map):
    """
    Simplify a K-expression.
    :return: a symbolic expression
    """
    assert len(expr) > 0
    if len(expr) == 1:  # must be a single terminal
        t = expr[0]
        assert isinstance(
            t, Terminal), 'A K-expression of length 1 must only contain a terminal.'
        if t.value is None:  # an input
            return sp.Symbol(t.name)
        return t.value

    expr = expr[:]  # because we need to change expr
    # K-expression is simply a level-order serialization of an expression tree.
    for i in reversed(range(len(expr))):
        p = expr[i]
        if isinstance(p, Function):
            try:
                sym_func = symbolic_function_map[p.name]
            except KeyError as e:
                print(
                    "Please provide the symbolic function mapping for '{}' in symbolic_function_map.".format(p.name))
                raise e
            args = []
            for _ in range(p.arity):
                t = expr.pop()  # t may be a terminal or a symbolic expression already
                if isinstance(t, Terminal):
                    if isinstance(t, SymbolTerminal):
                        args.append(sp.Symbol(t.name))
                    else:
                        args.append(t.value)
                else:
                    args.append(t)
            # evaluate this function node symbolically
            r = sym_func(*reversed(args))
            expr[i] = sp.simplify(r)
    return expr[0]


def simplify(genome, symbolic_function_map=None):
    """
    Compile the primitive tree into a (possibly simplified) symbolic expression.

    :param genome: :class:`~geppy.core.entity.KExpression`, :class:`~geppy.core.entity.Gene`, or
        :class:`~geppy.core.entity.Chromosome`, the genotype of an individual
    :param symbolic_function_map: dict, maps each function name in the primitive set to a symbolic version
    :return: a (simplified) symbol expression

    For example, *add(sub(3, 3), x)* may be simplified to *x*. This :func:`simplify` function can be used to
    postprocess the best individual obtained in GEP for a simplified representation. Some Python functions like
    :func:`operator.add` can be used directly in *sympy*. However, there are also functions that have their own
    symbolic versions to be used in *sympy*, like the :func:`operator.and_`, which should be replaced by
    :func:`sympy.And`. In such a case, we may provide a map
    ``symbolic_function_map={operator.and_.__name__, sympy.And}`` supposing the function primitive encapsulating
    :func:`operator.and_` uses its default name.

    Such simplification doesn't affect GEP at all. It should be used as a postprocessing step to simplify the final
    solution evolved by GEP.

    .. note::
        If the *symbolic_function_map* argument remains as the default value ``None``, then a default map
        :data:`DEFAULT_SYMBOLIC_FUNCTION_MAP` is used, which contains common
        *name-to-symbolic function* mappings, including the arithmetic operators and Boolean logic operators..

    .. note::
        This function depends on the :mod:`sympy` module. You can find it `here <http://www.sympy.org/en/index.html>`_.
    """
    if symbolic_function_map is None:
        symbolic_function_map = DEFAULT_SYMBOLIC_FUNCTION_MAP
    if isinstance(genome, KExpression):
        return _simplify_kexpression(genome, symbolic_function_map)
    elif isinstance(genome, Gene):
        return _simplify_kexpression(genome.kexpression, symbolic_function_map)
    elif isinstance(genome, Chromosome):
        if len(genome) == 1:
            return _simplify_kexpression(genome[0].kexpression, symbolic_function_map)
        else:   # multigenic chromosome
            simplified_exprs = [_simplify_kexpression(
                g.kexpression, symbolic_function_map) for g in genome]
            # combine these sub-expressions into a single one with the linking function
            try:
                linker = symbolic_function_map[genome.linker.__name__]
            except:
                linker = genome.linker
            return sp.simplify(linker(*simplified_exprs))
    else:
        raise TypeError('Only an argument of type KExpression, Gene, and Chromosome is acceptable. The provided '
                        'genome type is {}.'.format(type(genome)))
