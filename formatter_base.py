from __future__ import annotations
from abc import ABC
from functools import singledispatchmethod
from chain import Function


class Formatter(ABC):
    """Abstract formatter. Subclasses should register handlers for
    concrete Function subclasses using Formatter.format.register(Type).

    The default implementation raises NotImplementedError so missing
    handlers are easy to detect.
    """

    @singledispatchmethod
    def format(self, f: Function) -> str:  # default for unknown types
        raise NotImplementedError(f"No formatter registered for type {type(f)}")
