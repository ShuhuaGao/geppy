# coding=utf-8
"""
The *geppy* package supports gene expression programming in Python.

.. tip::
    For convenience of usage, all APIs have been imported into the top-level package *geppy* and can be used with the
    short form ``geppy.API`` directly. For example, simply use ``geppy.Chromosome`` as an alias of the class
    :class:`geppy.core.entity.Chromosome` and ``geppy.invert`` as the alias of the function
    :func:`geppy.tools.mutation.invert`.
"""
#    This file is part of geppy.
#
#    geppy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    geppy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with geppy. If not, see <http://www.gnu.org/licenses/>.

from pkg_resources import get_distribution, DistributionNotFound

# fetch version from setup.py
try:
    __version__ = get_distribution('geppy').version
except DistributionNotFound as e:
    __version__ = 'Please install this package with setup.py'

__author__ = 'Shuhua Gao'
# __revision__ = "1.2.2"

from .core.entity import *
from .core.symbol import *
from .tools.crossover import *
from .tools.mutation import *
from .algorithms.basic import *
from .tools.parser import *
from .tools.generator import *
from .support.visualization import *
from .tools.toolbox import *


def _print_module_not_found(module_name, f_name, http=''):
    print('''The module {0} is not available, on which {1} depends.
                The geepy.graph function will not be visible until this module is installed. 
                To install graphviz, please check {2}.'''.format(module_name, f_name, http))


try:
    import graphviz
    from geppy.support.visualization import export_expression_tree
except ImportError:
    _print_module_not_found('graphviz', 'geepy.export_expression_tree', 'https://pypi.org/project/graphviz/')

try:
    import sympy
    from geppy.support.simplification import simplify
except ImportError:
    _print_module_not_found('sympy', 'geepy.simplify', 'http://www.sympy.org/')


