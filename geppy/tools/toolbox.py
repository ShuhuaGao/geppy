# coding=utf-8
"""
.. moduleauthor:: Shuhua Gao

The module :mod:`toolbox` mainly provides a class :class:`Toolbox` on the basis of :class:`deap.base.Toolbox` for GEP
operator registration, which is then extensively used in the builtin GEP algorithms. For user-defined custom algorithms,
they should also accept a toolbox instance and retrieve the operators and their associated probability from the toolbox
instance.
"""

from deap import base
import collections


class Toolbox(base.Toolbox):
    """
    A toolbox for evolution that contains the evolutionary operators.
    Initially, the toolbox contains a :meth:`~Toolbox.clone` method that
    duplicates any element it is passed as argument, this method defaults to
    the :func:`copy.deepcopy` function. and a :meth:`map`
    method that applies the function given as first argument to every items
    of the iterables given as next arguments, this method defaults to the
    :func:`map` function. You may populate the toolbox with any other
    function by using the :meth:`~Toolbox.register` method.

    As an extension of :class:`deap.base.Toolbox`, this class adds a :attr:`~Toolbox.pbs` property to specify the
    probabilities for the registered operators. For example, `pbs['mut_uniform']` gives the probability for an
    operator with the alias `'mut_uniform'`. Besides, :attr:`Toolbox.pbs` is an :class:`~collections.OrderedDict` instance
    which can remember the order of entry insertion. Thus, if `'mut_A'` is inserted into :attr:`~Toolbox.pbs` before
    `'mut_B'`, then the operator with the alias `'mut_A'` is applied earlier than the one corresponding to `'mut_B'`.

    A short way to combine operator registration and probability specification together is to pass a keyword-only argument
    `pb` into the :meth:`register` method, which will insert the probability automatically into :attr:`pbs`.
    """

    def __init__(self):
        """
        Initialize an *empty* toolbox with only the default :meth:`clone` and :meth:`map` methods.
        """
        super().__init__()
        self._pbs = collections.OrderedDict()

    def register(self, alias, function, *args, **kargs):
        """
        Register a *function* in the toolbox under the name *alias*. You
        may provide default arguments that will be passed automatically when
        calling the registered function. Fixed arguments can then be overriden
        at function call time.

        :param alias: The name the operator will take in the toolbox. If the
                              alias already exists it will overwrite the the operator
                              already present.
        :param function: The function to which the alias refers.
        :param args: one or more positional arguments to pass to the registered function, optional
        :param kargs: one or more keyword arguments to pass to the registered function, optional

        .. hint::
            Under the hood lies the partial function binding. Check :func:`functools.partial` for details.

        .. note::
            If an operator needs its probability specified, like mutation and crossover operators, it can be done by
            inserting the probability into the :attr:`pbs` dictionary with the same alias. Alternatively, it can be
            given with a the special keyword argument `pb` in this method ::

                tb = Toolbox()
                tb.register('mut_uniform', mutate_uniform, ind_pb=0.02)
                tb.pbs['mut_uniform'] = 0.1

            or equivalently ::

                tb = Toolbox()
                tb.register('mut_uniform', mutate_uniform, ind_pb=0.02, pb=0.1)

            As a result, the special keyword argument `pb` is always excluded from binding into *function*.
        """
        pb = None
        if 'pb' in kargs:
            pb = kargs['pb']
            del kargs['pb']
        super().register(alias, function, *args, **kargs)
        if pb is not None:
            self.pbs[alias] = pb

    @property
    def pbs(self):
        """
        Get the probability dictionary of type :class:`~collections.OrderedDict`,
        which specifies the probabilities for the registered operators.
        """
        return self._pbs


__all__ = ['Toolbox']
