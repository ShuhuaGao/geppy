# coding=utf-8
"""
.. moduleauthor:: Shuhua Gao

This module :mod:`visualization` provides utility functions to visualization the expression tree from a given
K-expression, a gene or a chromosome in GEP.
"""

from ..core.entity import KExpression, Chromosome, Gene
from ..core.symbol import Function, Terminal


def _graph_kexpression(expr, starting_index):
    """
    Create a graph for a K-expression *expr* with the node's number starting from *starting_index*.

    :param expr: k-expression
    :param starting_index: the first number of nodes in the expression tree
    :return: A node list, an edge list, and a dictionary of labels.
    """
    assert len(expr) > 0
    nodes = [starting_index + i for i in range(len(expr))]
    edges = []
    labels = {}
    for i, p in enumerate(expr):
        if isinstance(p, Function):
            labels[starting_index + i] = p.name
        elif isinstance(p, Terminal):
            labels[starting_index + i] = p.format()
        else:
            raise RuntimeError('Unrecognized symbol. Normally, a symbol in the K-expression is either a function '
                               'or a terminal')
    i = 0
    j = 0
    while i < len(expr):
        for _ in range(expr[i].arity):
            j += 1
            edges.append((i + starting_index, j + starting_index))
        i += 1
    return nodes, edges, labels


def graph(genome, label_renaming_map=None):
    """
    Construct the graph of a genome. It returns in order a node list, an edge list, and a dictionary of the per node
    labels. The node are represented by numbers, the edges are tuples connecting two nodes (number), and the labels are
    values of a dictionary for which keys are the node numbers.

    :param genome: :class:`~geppy.core.entity.KExpression`, :class:`~geppy.core.entity.Gene`, or
        :class:`~geppy.core.entity.Chromosome`, the genotype of an individual
    :param label_renaming_map: dict, which maps the old name of a primitive (or a linking function)
        to a new one for better visualization. The default label for each node is just the name of the primitive
        placed on this node. For example, you may provide ``renamed_labels={'and_': 'and'}``.
    :return: A node list, an edge list, and a dictionary of labels.

    You can visualize a genome and export the tree visualization to an image file directly using the
    :func:`export_expression_tree` function.
    """
    nodes = []
    edges = []
    labels = {}
    if isinstance(genome, KExpression):
        nodes, edges, labels = _graph_kexpression(genome, 0)
    elif isinstance(genome, Gene):
        nodes, edges, labels = _graph_kexpression(genome.kexpression, 0)
    elif isinstance(genome, Chromosome):
        if len(genome) == 1:
            nodes, edges, labels = _graph_kexpression(genome[0].kexpression, 0)
        else:   # multigenic chromosome, we need to concatenate multiple trees
            starting_index = 1
            sub_roots = []
            for gene in genome:
                expr = gene.kexpression
                sub_roots.append(starting_index)
                sub_nodes, sub_edges, sub_labels = _graph_kexpression(
                    expr, starting_index)
                nodes.extend(sub_nodes)
                edges.extend(sub_edges)
                labels.update(sub_labels)
                starting_index += len(expr)
            # connect subtrees by inserting the linker node as 0
            nodes.append(0)
            for root in sub_roots:
                edges.append((0, root))
            labels[0] = genome.linker.__name__
    else:
        raise TypeError('Only an argument of type KExpression, Gene, and Chromosome is acceptable. The provided '
                        'genome type is {}.'.format(type(genome)))
    # rename_labels labels
    if label_renaming_map is not None:
        for k, v in labels.items():
            if v in label_renaming_map:
                labels[k] = label_renaming_map[v]
    return nodes, edges, labels


def export_expression_tree(genome, label_renaming_map=None, file='tree.png'):
    """
    Construct the graph of a *genome* and then export it to a *file*.

    :param genome: :class:`~geppy.core.entity.KExpression`, :class:`~geppy.core.entity.Gene`, or
        :class:`~geppy.core.entity.Chromosome`, the genotype of an individual
    :param label_renaming_map: dict, which maps the old name of a primitive (or a linking function)
        to a new one for better visualization. The default label for each node is just the name of the primitive
        placed on this node. For example, you may provide ``renamed_labels={'and_': 'and'}``.
    :param file: str, the file path to draw the expression tree, which may be a relative or absolute one.
        If no extension is included in *file*, then the default extension 'png' is used.

    .. note::
        This function currently depends on the :mod:`graphviz` module to render the tree. Please first install the
        `graphviz <https://pypi.org/project/graphviz/>`_ module before using this function.
        Alternatively, you can always obtain the raw graph data with the :func:`graph` function, then postprocess the
        data and render them with other tools as you want.
    """
    import graphviz as gv
    import os.path

    _, edges, labels = graph(genome, label_renaming_map)
    file_name, ext = os.path.splitext(file)
    ext = ext.lstrip('.')
    g = gv.Graph(format=ext)
    for name, label in labels.items():
        g.node(str(name), str(label))  # add node
    for name1, name2 in edges:
        g.edge(str(name1), str(name2))  # add edge
    g.render(file_name)


__all__ = ['graph', 'export_expression_tree']
