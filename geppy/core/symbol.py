# coding=utf-8
"""
.. moduleauthor:: Shuhua Gao

This :mod:`symbol` module defines the classes encapsulating the symbols in GEP. A gene in GEP is composed of multiple
symbols,  which include *terminals* and *functions*. Together, they are called *primitives*. This module also
provides a :class:`PrimitiveSet` class, which represents a primitive set containing :class:`Terminal`
and :class:`Function` items, for management purpose.

This module implementation refers a lot to the
`deap.gp <http://deap.readthedocs.io/en/master/tutorials/advanced/gp.html>`_ module in DEAP,
especially on primitive management and parsing during evaluation. However, unlike genetic programming (GP) in
:mod:`deap.gp`, in GEP we usually only consider loosely  typed evolution, i.e., the terminals and functions don't
require specific data types. As a result, the codes can be greatly simplified.

.. note::
    In GEP and GP terminology, there are two kinds of primitives: *functions* and *terminals*. In
    `deap.gp <http://deap.readthedocs.io/en/master/tutorials/advanced/gp.html>`_, a *primitive* actually refers to
    a function. To avoid possible ambiguity, in :mod:`geppy` they are explicitly named *function* and *terminal*
    respectively. For instance, the :meth:`PrimitiveSet.add_function` method can be used to add a function primitive.

**Reference:**

.. [FC2006] Ferreira, Cândida. Gene expression programming: mathematical modeling by an artificial
                intelligence. Vol. 21. Springer, 2006.
"""
import numbers
import keyword


def _is_identifier_not_keyword(var):
    """
    Check whether the string *var* is a valid non-keyword Python identifier.
    For example, 'ab' and 'x1' satisfy this condition, while 'if' and 'lambda' are keywords.
    :param var: a string
    :return: bool
    """
    if var.isidentifier() and not keyword.iskeyword(var):
        return True
    return False


class Primitive:
    """
    Class that encapsulates a primitive in GEP. A primitive may be a function or a terminal.
    """
    def __init__(self, name, arity):
        """
        Initialize a primitive.

        :param name: str, name of the primitive
        :param arity: int, arity of the primitive
        """
        self._name = name
        self._arity = arity

    @property
    def arity(self):
        """
        Get the arity of this primitive. For a terminal, the arity is 0.
        """
        return self._arity

    @property
    def name(self):
        """
        Get the name of this primitive.
        """
        return self._name

    def format(self, *args):
        """
        Get a string representation of the underlying concept with possible arguments. The derived classes should
        implement this method.

        :param args: variable number of arguments
        :return: str, a string
        """
        raise NotImplementedError()

    def __str__(self):
        return self.name

    def __repr__(self):
        return '{}(name={}, arity={})'.format(self.__class__, self.name, self.arity)


class Function(Primitive):
    """
    Class that encapsulates a function in GEP. Note that this class only stores the function ID, i.e., the *name*
    attributed instead of the callable function itself. On the tmp hand, the underlying callable is retrieved
    somewhere else when needed, for example, from :class:`PrimitiveSet`. Thus, in the whole GEP program the provided
    name for each function must be unique.
    """
    def __init__(self, name, arity):
        """
                Initialize a function.

                :param name: str, name of the function
                :param arity: int, arity of the function
                """
        super().__init__(name, arity)
        args = ', '.join('{{{}}}'.format(index) for index in range(arity))  # like '{0}, {1}, {2}'
        self._seq = name + '(' + args + ')'  # e.g., add, --> 'add({0}, {1})'

    def format(self, *args):
        """
        Insert the arguments *args* into the function and get a Python statement to call the functions with the
        arguments in a string form. This returned string can afterwards by evaluated using the builtin
        :func:`eval` function.

        >>> f = Function('add', 2)
        >>> f.format(2, 10)
        'add(2, 10)'

        :param args: arguments, whose number should be equal to the arity of this function
        :return: str, a string form of function calling with arguments
        """
        assert len(args) == self.arity, "Function {} requires {} arguments while {} are provided.".format(
            self.name, self.arity, len(*args))
        return self._seq.format(*args)


class Terminal(Primitive):
    """
    Class that encapsulates a terminal in GEP.
    """
    def __init__(self, name, value):
        """
        Initialize a terminal.

        :param name: str, name of the terminal
        :param value: value of the terminal
        """
        super().__init__(name, arity=0)
        self._value = value

    def format(self):
        """
        Get a string representation of this terminal for subsequent evaluation purpose, which may be a symbol or a
        numeric string.

        :return: str, a string representation
        """
        return self.name

    @property
    def value(self):
        """
        Get the value of this terminal.
        """
        return self._value

    def __repr__(self):
        return '{}(name={}, value={})'.format(self.__class__, self.name, self.value)


class Ephemeral(Terminal):
    """
    Class that encapsulates an ephemeral random constant terminal in GEP, which can be used to build a normal constant
    terminal with a randomly produced value when needed using the :meth:`Ephemeral.generate` method.

    .. note::
            This special terminal named *ephemeral random constant* was originally introduced in genetic programming
            by Koza to handle numerical constants in evolutionary programming. However, it is NOT the recommended way
            to evolve constants in GEP, especially for complex problems. The recommended method is to add a RNC domain
            in genes.See Chapter 5 of [FC2006]_.
    """
    def __init__(self, name, gen):
        """

        :param name: str, name of the ephemeral constant
        :param gen: callable, an ephemeral number generator, which should return a random value when called with
            no arguments.
        """
        super().__init__(name, value=None)
        self._gen = gen

    @property
    def generator(self):
        """
        Get the ephemeral random number generator.
        """
        return self._gen

    @property
    def value(self):
        """
        This doesn't make sense for an ephemeral constant terminal. Simply raise :class:`NotImplementedError`.
        """
        raise NotImplementedError

    @property
    def format(self):
        """
        This doesn't make sense for an ephemeral constant terminal. Simply raise :class:`NotImplementedError`.
        """
        raise NotImplementedError

    def generate(self):
        """
        Generate a constant terminal, whose value is randomly produced by the :meth:`~Ephemeral.generator`.
        """
        value = self._gen()
        return Terminal(str(value), value)

    def __repr__(self):
        return '{}(name={}, gen={})'.format(self.__class__, self.name, self.generator)


class TerminalRNC(Terminal):
    """
    A special terminal, which is just a placeholder representing a random numerical constant (RNC) in the GEP-RNC
    algorithm. This class is mainly used internally by the :class:`~geppy.core.entity.GeneDc` class. The name of a
    :class:`TerminalRNC` object is '?' by default, and its value is retrieved dynamically according to the GEP-RNC
    algorithm. Refer to Chapter 5 of [FC2006]_ for more details.
    """
    def __init__(self, name='?', value=None):
        """
        Initialize a RNC terminal.

        :param name: str, name of the terminal
        :param value: value of the terminal
        """
        super().__init__(name, value)


class PrimitiveSet:
    """
    A class representing a primitive set, which contains the primitives (terminals and functions) that are used in GEP.

    .. warning::
        Each primitive in this set must have a unique name, which should also be a valid non-keyword
        Python identifier, except the constant terminals and RNC terminals. For more details, refer to the
        :meth:`~PrimitiveSet.add_terminal` and :meth:`add_rnc` methods.

    .. note::
        To use the GEP-RNC algorithm, i.e., to use :class:`~geppy.core.entity.GeneDc` for numerical constant handling,
        please use the :meth:`PrimitiveSet.add_rnc` method, which will add a special terminal
        of type :class:`~geppy.core.symbol.TerminalRNC` internally. Then, use a chromosome composed of
        :class:`~geppy.core.entity.GeneDc` genes as the individual.
        Refer to Chapter 5 of [FC2006]_ for more details.
    """
    def __init__(self, name, input_names):
        """
        Initiate a primitive set with the given *name* and the list of input names.

        :param name: name of the primitive set
        :param input_names: iterable, a list of names for the inputs in the GEP problem, for instance,
            ``['x', 'y']`` for two inputs.
        """
        self._name = name
        self._functions = []
        self._terminals = []
        self._globals = {'__builtins__': None}
        self._inputs = input_names
        self._unique_names = set()
        # add input terminals
        for name in input_names:
            self._add_input_terminal(name)

    def add_function(self, func, arity, name=None):
        """
        Add a function, which is internally encapsulated as a :class:`Function` object.

        :param func: callable
        :param arity: number of arguments accepted by *func*
        :param name: name of *func*. If remaining ``None``, then the ``func.__name__`` attribute is used.
        """
        if name is None:
            name = func.__name__
        self._assert_symbol_valid_and_unique(name)
        function_ = Function(name, arity)
        self._functions.append(function_)
        self._globals[name] = func

    def _add_input_terminal(self, name):
        """
        Add an input terminal with the given *name*.
        :param name: str, which must be a valid Python identifier.
        """
        self._assert_symbol_valid_and_unique(name)
        t = Terminal(name, value=None)
        self._terminals.append(t)

    def _add_constant_terminal(self, value):
        """
        Add a terminal which is a constant. There is no need to store such kind of terminals into *globals* for later
        evaluation, because their string representation acquired by :meth:`~Terminal.format` can be used directly.
        :param value: value of the terminal. Only numeric and Boolean types can be accepted.
        """
        assert isinstance(value, numbers.Number) or isinstance(value, bool),\
            '''Only a number or a Boolean value can be used to create a constant terminal. The provided value is {} of 
                type {}'''.format(value, type(value))
        self._terminals.append(Terminal(str(value), value))

    def _add_symbolic_terminal(self, value, name):
        """
        Add a symbolic terminal, whose name differs from the value. For example, we may add a terminal for the constant
        pi, with ``name='pi'`` and ``value=3.1415926``. This kind of terminals are also inserted into *globals* for
        later evaluation.

        :param value: value of the terminal
        :param name: name fo the terminal
        """
        self._assert_symbol_valid_and_unique(name)
        t = Terminal(name, value)
        self._globals[name] = value
        self._terminals.append(t)

    def _assert_symbol_valid_and_unique(self, name):
        """
        Assert whether a *name* is a valid symbol in this set. It should be (1) unique in *globals* and (2) a valid
        Python identifier. If the supplied *name* is acceptable, it is added to the :attr:`~PrimitiveSet._unique_names`
        attribute.
        """
        assert _is_identifier_not_keyword(name), 'A name must be a valid non-keyword Python identifier.' \
                                                  "The provided name '{}' cannot be accepted.".format(name)
        assert name not in self._unique_names, "A primitive with the name '{}' already exists.".format(name)
        self._unique_names.add(name)

    def add_terminal(self, value, name=None):
        """
        Add a non-ephemeral terminal, which is internally encapsulated as a :class:`Terminal` object.

        :param value: value of the terminal
        :param name: str, name of the terminal, default: ``None``.

        If a constant terminal like a number 5 is expected, then simply leave *name* as ``None``, and use
        ``add_terminal(5)``. Otherwise, please specify a valid *name* to create a symbolic terminal.
        For instance, ``add_terminal(3.14, 'pi')``.
        """
        if name is None:
            self._add_constant_terminal(value)
        else:
            self._add_symbolic_terminal(value, name)

    def add_ephemeral(self, gen, name):
        """
        Add an ephemeral constant to the set. An ephemeral's true value is generated by *gen*, which should be a no
        argument function that returns a random value, when the ephemeral is picked to build a gene. Once an ephemeral
        is picked, its :meth:`Ephemeral.generate` method is called to produce a constant terminal with a random value
        generated by *gen*. Thus, each time :meth:`Ephemeral.generate` is called, a terminal with a different value may
        be yielded.

        :param gen: callable, a random number generator which returns a random constant when called by ``gen()``
        :param name: str, name of the ephemeral

        .. note::
            The special terminal named *ephemeral random constant* was originally introduced in genetic programming
            by Koza to handle numerical constants in evolutionary programming. However, though an ephemeral still works
            in GEP, it is NOT the recommended way to evolve constants in GEP, especially for complex problems.
            The recommended method is to add a RNC domain in genes. See Chapter 5 of [FC2006]_.
        """
        self._assert_symbol_valid_and_unique(name)
        self._terminals.append(Ephemeral(name, gen))

    def add_rnc(self, name='?'):
        """
        Add a special terminal representing a random numerical constant (RNC), as defined in the GEP-RNC algorithm.
        This terminal's value is retrieved dynamically from an RNC array attached to a gene of type
        :class:`~geppy.core.entity.GeneDc` according to the GEP-RNC algorithm.
        See also :class:`TerminalRNC` and refer to Chapter 5 of [FC2006]_ about GEP-RNC.

        :param name: str, name of the terminal. For a RNC terminal, generally there is no need to specify a name and
            the default value '?' is recommended.

        Usually it is sufficient to call this method once to add only one RNC terminal, since the values of the
        RNC terminals at different positions are different random numerical constants..
        """
        self._terminals.append(TerminalRNC(name))

    @property
    def functions(self):
        """
        Get all the functions.
        """
        return self._functions

    @property
    def terminals(self):
        """
        Get all the terminals, including the ephemeral ones and input terminals if any.
        """
        return self._terminals

    @property
    def ephemerals(self):
        """
        Get all the ephemeral terminals.
        """
        return [t for t in self.terminals if isinstance(t, Ephemeral)]

    @property
    def name(self):
        """
        Get the name of this primitive set.
        """
        return self._name

    @property
    def input_names(self):
        """
        Get a list of names for the input arguments.
        """
        return self._inputs

    def __str__(self):
        terminals = '[{}]'.format(', '.join(str(t) for t in self.terminals))
        functions = '[{}]'.format(', '.join(str(f) for f in self.functions))
        return 'PrimitiveSet {}\n\tFunctions: {}\n\tTerminals: {}'.format(self.name, functions, terminals)

    @property
    def globals(self):
        """
        Get a dictionary which can be used to set up the evaluation/execution environment. This dictionary can be
        fed into the builtin :func:`eval` and :func:`exec` functions for expression evaluation.

        For example, if we call ``add_function(max, 2, 'max2')``, then ``globals['max2']`` corresponds to a
        :class:`Function` object encapsulating the :func:`max` function.
        """
        return self._globals

    @property
    def max_arity(self):
        """
        Get the max arity of functions in this primitive set.
        """
        return max(f.arity for f in self.functions)
