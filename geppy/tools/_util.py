"""
Some utility functions used internally.
"""

import random
import copy
from ..core.symbol import EphemeralTerminal


def _choose_a_terminal(terminals):
    """
    Choose a terminal from the given list *terminals* randomly.

    :param terminals: iterable, terminal set
    :return: a terminal
    """
    terminal = random.choice(terminals)
    if isinstance(terminal, EphemeralTerminal):  # an Ephemeral
        terminal = copy.deepcopy(terminal)  # generate a new one
        terminal.update_value()
    return terminal
