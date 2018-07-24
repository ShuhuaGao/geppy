# coding=utf-8
"""
.. module:: geppy.core.entity
.. moduleauthor:: Shuhua Gao

The module :mod:`entity` defines the core data structures used in GEP, including the gene, the chromosome, the
K-expression and the expression tree. We use primitives (functions and terminals) to build genes and then compose
a chromosome with one or multiple genes. Refer to the :mod:`geppy.core.symbol` module for information on primitives.

A chromosome composed of only one gene is called a *monogenic* chromosome, while one composed of more than one genes
is named a *multigenic* chromosome. A multigenic chromosome can be assigned a linking function to combine the results
from multiple genes into a single result.
"""

from ..tools.generator import generate_genome
from ..core.symbol import Function


_DEBUG = False


class KExpression(list):
    """
    Class representing the K-expression, or the open reading frame (ORF), of a gene in GEP. A K-expression is actually
    a linear form of an expression tree obtained by level-order traversal.
    """
    def __init__(self, content):
        """
        Initialize a K-expression.

        :param content: iterable, each element is a primitive (function or terminal)
        """
        list.__init__(self, content)

    def __str__(self):
        return ', '.join(ele.name for ele in self)


class ExpressionTree:
    """
    Class representing an expression tree (ET) in GEP, which may be obtained by translating a K-expression, a
    gene, or a chromosome, i.e., genotype-phenotype mapping.
    """
    def __init__(self, root):
        """
        Initialize a tree with the given *root* node.

        :param root: :class:`ExpressionTree.Node`, the root node
        """
        self._root = root

    @property
    def root(self):
        """
        Get the root node of this expression tree.
        """
        return self._root

    class Node:
        """
        Class representing a node in the expression tree. Each node has a variable number of children, depending on
        the arity of the primitive at this node.
        """
        def __init__(self, name):
            self._children = []
            self._name = name

        @property
        def children(self):
            """
            Get the children of this node.
            """
            return self._children

        @property
        def name(self):
            """
            Get the name (label) of this node.
            """
            return self._name

    @classmethod
    def from_genotype(cls, genome):
        """
        Create an expression tree by translating *genome*, which may be a K-expression, a gene, or a chromosome.

        :param genome: :class:`KExpression`, :class:`Gene`, or :class:`Chromosome`, the genotype of an individual
        :return: :class:`ExpressionTree`, an expression tree
        """
        if isinstance(genome, Gene):
            return cls._from_kexpression(genome.kexpression)
        elif isinstance(genome, KExpression):
            return cls._from_kexpression(genome)
        elif isinstance(genome, Chromosome):
            if len(genome) == 1:
                return cls._from_kexpression(genome[0].kexpression)
            sub_trees = [cls._from_kexpression(gene.kexpression) for gene in genome]
            # combine the sub_trees with the linking function
            root = cls.Node(genome.linker.__name__)
            root.children[:] = sub_trees
            return cls(root)
        raise TypeError('Only an argument of type KExpression, Gene, and Chromosome is acceptable. The provided '
                        'genome type is {}.'.format(type(genome)))

    @classmethod
    def _from_kexpression(cls, expr):
        """
        Create an expression tree from a K-expression.

        :param expr: a K-expression
        :return: :class:`ExpressionTree`, an expression tree
        """
        if len(expr) == 0:
            return None
        # first build a node for each primitive
        nodes = [cls.Node(p.name) for p in expr]
        # connect each node to its children if any
        i = 0
        j = 0
        while i < len(nodes):
            for _ in range(expr[i].arity):
                j += 1
                nodes[i].children.append(nodes[j])
            i += 1
        return cls(nodes[0])


class Gene(list):
    """
    A single gene in GEP, which is a fixed-length linear sequence composed of functions and terminals.
    """
    def __init__(self, pset, head_length):
        """
        Instantiate a gene.
        :param head_length: length of the head domain
        :param pset: a primitive set including functions and terminals for genome construction.

        Supposing the maximum arity of functions in *pset* is *max_arity*, then the tail length is automatically
        determined to be ``tail_length = head_length * (max_arity - 1) + 1``. The genome, i.e., list of symbols in
        the instantiated gene is formed randomly from *pset*.
        """
        self._head_length = head_length
        genome = generate_genome(pset, head_length)
        list.__init__(self, genome)

    @classmethod
    def from_genome(cls, genome, head_length):
        """
        Build a gene directly from the given *genome*.

        :param genome: iterable, a list of symbols representing functions and terminals
        :param head_length: length of the head domain
        :return: :class:`Gene`, a gene
        """
        g = super().__new__(cls)
        super().__init__(g, genome)
        g._head_length = head_length
        return g

    def __deepcopy__(self, memodict):
        """
        Deep copy a gene. Note that there is no need to deep copy the contained primitives since they will not be
        changed during GEP.
        If we rely on the default deepcopy implementation, then the contents are also deep copied unnecessarily.
        """
        return self.__class__.from_genome(self, head_length=self.head_length)

    def __str__(self):
        """
        Return the expression in a human readable string.

        :return: string form of the expression
        """
        expr = self.kexpression
        i = len(expr) - 1
        while i >= 0:
            if expr[i].arity > 0:  # a function
                f = expr[i]
                args = []
                for _ in range(f.arity):
                    ele = expr.pop()
                    if isinstance(ele, str):
                        args.append(ele)
                    else:  # a terminal or an ephemeral class
                        args.append(ele.format())
                expr[i] = f.format(*reversed(args))  # replace the operator with its result (str)
            i -= 1

        return expr[0] if isinstance(expr[0], str) else expr[0].format()

    @property
    def kexpression(self):
        """Get the K-expression represented by this gene.
        """
        # get the level-order K expression
        expr = KExpression([self[0]])
        i = 0
        j = 1
        while i < len(expr):
            for _ in range(expr[i].arity):
                expr.append(self[j])
                j += 1
            i += 1
        return expr

    @property
    def head_length(self):
        """
        Get the length of the head domain.
        """
        return self._head_length

    @property
    def tail_length(self):
        """
        Get the length of the tail domain.
        """
        return len(self) - self.head_length

    @property
    def max_arity(self):
        """
        Get the max arity of the functions in the primitive set with which the gene is built.
        """
        return (len(self) - 1) // self.head_length

    @property
    def head(self):
        """
        Get the head domain of the gene.
        """
        return self[0: self.head_length]

    @property
    def tail(self):
        """
        Get the tail domain of the gene.
        """
        return self[self.head_length:]

    def __repr__(self):
        """
        Return a string representation of a gene including all the symbols. Mainly for debugging purpose.

        :return: a string
        """
        info = []
        for p in self:
            if isinstance(p, Function):
                info.append(p.name)
            else:   # normally, a terminal or an ephemeral constant
                try:
                    info.append(p.format())
                except:
                    info.append(str(p))
        return '{} ['.format(self.__class__) + ', '.join(info) + ']'


class Chromosome(list):
    """
    Individual representation in gene expression programming (GEP), which may contain one or more genes. Note that in a
    multigenic chromosome, all genes must share the same head length and the same primitive set. Each element of
    :class:`Chromosome` is of type :class:`Gene`.
    """

    def __init__(self, gene_gen, n_genes, linker=None):
        """
        Initialize an individual in GEP.

        :param gene_gen: callable, a gene generator, i.e., ``gene_gen()`` should return a gene
        :param n_genes: number of genes
        :param linker: callable, a linking function

        .. note:
            If a linker is specified, then it must accept *n_genes* arguments, each produced by one gene. If the
            *linker* parameter is the default value ``None``, then for a monogenic chromosome, no linking is applied,
            while for a multigenic chromosome, the ``tuple`` function is used.
        """
        list.__init__(self, (gene_gen() for _ in range(n_genes)))
        self._linker = linker

    @classmethod
    def from_genes(cls, genes, linker=None):
        """
        Create a chromosome instance by providing a list of genes directly. The number of genes is implicitly specified
        by ``len(genes)``.

        :param genes: iterable, a list of genes
        :param linker: callable, a linking function. Refer to :meth:`Chromosome.__init__` for more details.
        :return: :class:`Chromosome`, a chromosome

        .. caution:
            Every gene of type :class:`Gene` in *genes* should be an independent object. That is, for any two genes
            *g1* and *g2*, the assertion ``g1 is g1`` should return ``False``. Otherwise, unexpected results may happen
            in following evolution.

        """
        c = super().__new__(cls)
        super().__init__(c, genes)
        c.linker = linker
        return c

    @property
    def head_length(self):
        """
        Get the length of the head domain. All genes in a chromosome should have the same head length.
        """
        return self[0].head_length

    @property
    def tail_length(self):
        """
        Get the length of the tail domain. All genes in a chromosome should have the same tail length.
        """
        return self[0].tail_length

    @property
    def max_arity(self):
        """
        Get the max arity of the functions in the primitive set with which all genes are built. Note that all genes
        in a chromosome should share the same primitive set.
        """
        return self[0].max_arity

    def __str__(self):
        """
        Return the expressions in a human readable string.
        """
        return '[\n\t' + ',\n\t'.join(str(g) for g in self) + '\n]'

    @property
    def kexpressions(self):
        """
        Get a list of K-expressions for all the genes in this chromosome.
        """
        return [g.kexpression for g in self]

    @property
    def linker(self):
        """
        Get the linking function.
        """
        if self._linker is None and len(self) > 1:
            return tuple
        return self._linker

    def __repr__(self):
        """
        Return a string representation of the chromosome including all the symbols in all the genes.
        Mainly for debugging purpose.

        :return: a string
        """
        return '{}[\n\t'.format(self.__class__) + ',\n\t'.join(repr(g) for g in self) + \
               '\n], linker={}'.format(self.linker)