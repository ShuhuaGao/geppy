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
import keyword
import copy
import numbers


def _is_nonkeyword_identifier(var):
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

    def __deepcopy__(self, memodict):
        """
        Faster deep copy because generally a primitive is immutable except the ephemeral terminals
        The derived classes should rewrite this method if needed.
        """
        return self


class Function(Primitive):
    """
    Class that encapsulates a function in GEP. Note that this class only stores the function ID, i.e., the *name*
    attribute instead of the callable function itself. On the other hand, the underlying callable is retrieved
    somewhere else when needed, for example, from :class:`PrimitiveSet`. Thus, in the whole GEP program the provided
    name for each function must be unique.
    """

    def __init__(self, name, arity):
        """
        Initialize a function.

        :param name: str, name of the function, must be valid non-keyword Python identifier
        :param arity: int, arity of the function
        """
        assert _is_nonkeyword_identifier(name), \
            'Name of a function must be a valid non-keyword Python identifier'
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

    def __str__(self):
        return self.format()


class ConstantTerminal(Terminal):
    """
    Class that represents a constant terminal, whose value will never change during the whole evolution.
    The value of a constant terminal can only be a literal variable , such as a numeric number or a Boolean variable.
    """
    def __init__(self, value):
        """
        Initialize a constant terminal with the given *value*.

        :param value: int, float, bool, value of the constant, can only be a number or a Boolean variable
        """
        assert isinstance(value, numbers.Number) or isinstance(value, bool), \
            'The value of a constant terminal can only be a number or a bool variable.'
        super().__init__(repr(value), value)


class SymbolTerminal(Terminal):
    """
    Class that represents a symbolic terminal. Only the :attr:`name` of a symbol terminal participates in the genome
    construction, while its value is retrieved from the :class:`PrimitiveSet` only during evaluation. Therefore, the
    :attr:`name` of a symbol terminal must be a valid and unique Python identifier.
    The value of a symbol terminal can be any object that leads to a legal expression when evaluated.
    """
    def __init__(self, name):
        """
        Initialize a symbol terminal with the symbol given by *name*.

        :param name: str, must be a valid non-keyword Python identifier
        """
        assert _is_nonkeyword_identifier(name), \
            'Name of a symbol terminal must be a valid non-keyword Python identifier'
        super().__init__(name, value='symbol')


class EphemeralTerminal(Terminal):
    """
    Class that encapsulates an ephemeral numeric constant terminal in GEP, which can be used to build a normal constant
    terminal with a randomly produced value when needed using the :meth:`Ephemeral.generate` method.

    Just as the name implies, the value of an ephemeral terminal may be changed during evolution, either by mutation or
    by some local search heuristics, to optimize the numerical constants in a mathematical model.

    .. note::
            This special terminal named *ephemeral random constant* was originally introduced in genetic programming
            by Koza to handle numerical constants in evolutionary programming. Another way to evolve simple constants in
            GEP is to add a RNC domain in genes and use the GEP-RNC algorithm. See Chapter 5 of [FC2006]_.

    .. seealso::
        :func:`~geppy.tools.mutation.mutate_uniform_ephemeral` for ephemeral mutation, and
        :class:`~geppy.core.entity.GeneDc` for the GEP-RNC algorithm.
    """

    def __init__(self, name, gen):
        """
        Initialize an ephemeral terminal with the given *name* and the random number generator *gen*. Its initial value
        is randomly generated by *gen* and can later be set explicitly with the :attr:`value` property. Alternatively,
        its value can be updated to a new value generated by *gen* with the :meth:`update_value` method. Note that just
        like a constant terminal, the value of an ephemeral terminal can only be a number or a bool.

        :param name: str, name of the ephemeral constant
        :param gen: callable, an ephemeral number generator, which should return a random value when called with
            no arguments ``gen()``.
        """
        super().__init__(name, value=None)
        self._gen = gen
        self.value = self._gen()

    @property
    def generator(self):
        """
        Get the ephemeral random number generator.
        """
        return self._gen

    @property
    def value(self):
        """
        Get the value.
        """
        return self._value

    @value.setter
    def value(self, val):
        """
        Set the value to be *val*.
        """
        assert isinstance(val, numbers.Number) or isinstance(val, bool), \
            'The value of a constant terminal can only be a number or a bool variable.'
        self._value = val

    def format(self):
        """
        Get a string representation of this terminal for subsequent evaluation purpose, which is equivalently
        ``repr(self.value)``.
        """
        return repr(self.value)

    def update_value(self):
        """
        Update the value of this ephemeral constant in place, whose new value is randomly produced by the random number
        generator :attr:`generator`.
        """
        self._value = self._gen()

    def __repr__(self):
        return '{}(name={}, gen={})'.format(self.__class__, self.name, self.generator)

    def __deepcopy__(self, memodict={}):
        """
        Deep copy the ephemeral and returns a new ephemeral terminal with the same value.
        """
        return copy.copy(self)


class RNCTerminal(Terminal):
    """
    A special terminal, which is just a placeholder representing a random numerical constant (RNC) in the GEP-RNC
    algorithm. This class is mainly used internally by the :class:`~geppy.core.entity.GeneDc` class. The name of a
    :class:`RNCTerminal` object is '?' by default, and its value is retrieved dynamically according to the GEP-RNC
    algorithm. Refer to Chapter 5 of [FC2006]_ for more details.

    In *geppy* and the GEP-RNC algorithm, the RNC terminal is just a placeholder and the default name '?' is
    recommended.
    """

    def __init__(self, name='?'):
        """
        Initialize a RNC terminal.

        :param name: str, default '?', name of the terminal
        """
        super().__init__(name, value=None)


class PrimitiveSet:
    """
    A class representing a primitive set, which contains the primitives (terminals and functions) that are used in GEP.

    .. note::
        Each function :class:`Function` and symbol terminal :class:`SymbolTerminal`
        must have their unique names. This is because internally their true value (or
        the callable for a function) will be be stored in a dictionary :attr:`globals` inside :class:`PrimitiveSet` and
        later be retrieved when compiling a genome. Consequently, the name must be a valid non-keyword Python identifier.
        To learn more about the underlying mechanism, refer to the Python builtin function :func:`compile`.

    .. note::
        To use the GEP-RNC algorithm, i.e., to use :class:`~geppy.core.entity.GeneDc` for numerical constant handling,
        please use the :meth:`PrimitiveSet.add_rnc_terminal` method, which will add a special terminal
        of type :class:`~geppy.core.symbol.RNCTerminal` internally. Then, use a chromosome composed of
        :class:`~geppy.core.entity.GeneDc` genes as the individual.
        Refer to Chapter 5 of [FC2006]_ for more details.
    """

    def __init__(self, name, input_names):
        """
        Initiate a primitive set with the given *name* and the list of input names.

        :param name: name of the primitive set
        :param input_names: iterable, a list of names for the inputs in the GEP problem, for instance,
            ``['x', 'y']`` for two inputs. Internally, a symbol terminal is built for each input.
        """
        self._name = name
        self._functions = []
        self._terminals = []
        self._globals = {'__builtins__': None}
        self._inputs = input_names
        self._unique_names = set()
        # add input terminals
        for name in input_names:
            self._assert_name_unique(name)
            self.add_symbol_terminal(name, value='input')

    def add_function(self, func, arity, name=None):
        """
        Add a function, which is internally encapsulated as a :class:`Function` object.

        :param func: callable
        :param arity: number of arguments accepted by *func*
        :param name: name of *func*, default ``None``.
            If remaining ``None``, then the ``func.__name__`` attribute is used instead.
        """
        if name is None:
            name = func.__name__
        self._assert_name_unique(name)
        function_ = Function(name, arity)
        self._functions.append(function_)
        self._globals[name] = func

    def add_constant_terminal(self, value):
        """
        Add a terminal which is a constant.

        There is no need to store such kind of terminals into *globals* for later
        evaluation, because their string representation acquired by :meth:`~Terminal.format` can be used directly.

        :param value: value of the terminal. Only numeric and Boolean types can be accepted.
        """
        self._terminals.append(ConstantTerminal(value))

    def add_symbol_terminal(self, name, value):
        """
        Add a symbolic terminal, whose name points to the true value stored in :attr:`globals`.

        For example, we may add a terminal for the constant *pi*, with ``add_symbol_terminal(3.1415, 'pi')``.

        :param name: name of the terminal
        :param value: value of the terminal
        """
        self._assert_name_unique(name)
        t = SymbolTerminal(name)
        self._globals[name] = value
        self._terminals.append(t)

    def _assert_name_unique(self, name):
        """
        Assert whether a *name* is a valid symbol in this set. It should be (1) unique in *globals* and (2) a valid
        Python identifier. If the supplied *name* is acceptable, it is added to the :attr:`~PrimitiveSet._unique_names`
        attribute.
        """
        assert name not in self._globals, "A primitive with the name '{}' already exists.".format(name)

    def add_ephemeral_terminal(self, name, gen):
        """
        Add an ephemeral constant to the set. An ephemeral's true value is generated by *gen*, which should be a no
        argument function that returns a random value, when the ephemeral is picked to build a gene. Once an ephemeral
        is picked, its :meth:`Ephemeral.generate` method is called to produce a random value.

        In GEP, the ephemeral terminals (i.e., their values) can be evolved by simple mutation or advanced
        local search optimization.

        :param name: str, name of the ephemeral
        :param gen: callable, a random number generator which returns a random constant when called by ``gen()``

        Usually calling :meth:`add_ephemeral_terminal` once is sufficient to generate the numerical coefficients in the
        model, because during evolution multiple independent copies of the ephemeral terminal will be built with
        various values.
        Nevertheless, if the numerical coefficients have various ranges, it may be more efficient to add multiple
        ephemeral terminals.

        .. note::
            The special terminal named *ephemeral random constant* was originally introduced in genetic programming
            by Koza to handle numerical constants in evolutionary programming. Alternatively, to evolve the numerical
            constants in GEP, we can also use the GEP-RNC algorithm.
            See Chapter 5 of [FC2006]_.

        .. seealso::
            :func:`~geppy.tools.mutation.mutate_uniform_ephemeral` for ephemeral mutation, and
            :class:`~geppy.core.entity.GeneDc` for the GEP-RNC algorithm.
        """
        self._assert_name_unique(name)
        self._terminals.append(EphemeralTerminal(name, gen))

    def add_rnc_terminal(self, name='?'):
        """
        Add a special terminal representing a random numerical constant (RNC), as defined in the GEP-RNC algorithm.
        This terminal's value is retrieved dynamically from an RNC array attached to a gene of type
        :class:`~geppy.core.entity.GeneDc` according to the GEP-RNC algorithm.
        See also :class:`RNCTerminal` and refer to Chapter 5 of [FC2006]_ about GEP-RNC.

        :param name: str, name of the terminal. For a RNC terminal, generally there is no need to specify a name and
            the default value '?' is recommended.

        Usually it is sufficient to call this method once to add only one RNC terminal, since the values of the
        RNC terminals at different positions are different random numerical constants..
        """
        self._terminals.append(RNCTerminal(name))

    @property
    def functions(self):
        """
        Get all the functions.
        """
        return self._functions

    @property
    def terminals(self):
        """
        Get all the terminals.
        """
        return self._terminals

    @property
    def ephemerals(self):
        """
        Get all the ephemeral terminals in this set.
        """
        return [t for t in self.terminals if isinstance(t, EphemeralTerminal)]

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
        """
        Gets an overview of the functions and terminals in this primitive set.
        """
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
