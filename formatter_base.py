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

    # Optional hooks for function definition emission.
    def explicit_function_definition(self, f: Function) -> str:
        """Return a full definition for ExplicitFunction-like objects.
        Subclasses may override to produce a language-specific definition.
        """
        raise NotImplementedError()

    def opaque_function_definition(self, f: Function) -> str:
        """Return a full definition for OpaqueFunction-like objects.
        """
        raise NotImplementedError()
