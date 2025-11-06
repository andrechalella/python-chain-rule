from __future__ import annotations
from abc import ABC
from functools import singledispatchmethod
from chain import Function


class Formatter(ABC):
    @singledispatchmethod
    def format(self, f: Function) -> str:
        raise NotImplementedError(f"No formatter registered for type {type(f)}")
