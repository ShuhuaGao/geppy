"""
Some utility functions used internally.
"""

import random
from ..core.symbol import Ephemeral


def _choose_a_terminal(terminals):
    """
    Choose a terminal from the given list *terminals* randomly.

    :param terminals: iterable, terminal set
    :return: a terminal
    """
    terminal = random.choice(terminals)
    if isinstance(terminal, Ephemeral):  # an Ephemeral
        terminal = terminal.clone_with_update()  # generate a normal constant terminal with a random value
    return terminal
