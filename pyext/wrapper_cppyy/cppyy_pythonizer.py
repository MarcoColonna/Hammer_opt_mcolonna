#**** This file is a part of the HAMMER library
#**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
#**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
#**** Please note the MCnet academic guidelines; see GUIDELINES for details

from enum import IntEnum

from cppyy.gbl import Hammer  # type: ignore[import-not-found]
from cppyy._pythonization import add_pythonization  # type: ignore[import-not-found]


# C++ enum values are not exposed as Python integers by standalone cppyy, so they
# can't be used with chr(), list indices, or integer-expecting builtins. These
# IntEnums mirror the C++ definitions and make both backends behave consistently.

class RecordType(IntEnum):  # enum RecordType : char (IOTypes.hh)
    UNDEFINED = ord('u')
    HEADER = ord('b')
    EVENT = ord('e')
    HISTOGRAM = ord('h')
    RATE = ord('r')
    HISTOGRAM_DEFINITION = ord('d')

class WTerm(IntEnum):  # enum class WTerm : uint8_t (IndexTypes.hh)
    COMMON = 0
    NUMERATOR = 1
    DENOMINATOR = 2

class PAction(IntEnum):  # enum class PAction : uint8_t (Hammer.hh)
    ALL = 0
    WEIGHTS = 1
    HISTOGRAMS = 2


def pythonize_four_momentum(klass, name):
    def __str__(self):
        return f"FourMomentum({self.E()}, {self.px()}, {self.py()}, {self.pz()})"

    def __mul__(self, other):
        if isinstance(other, (float, int)):
            res = Hammer.FourMomentum(self)
            res *= other
            return res
        elif isinstance(other, type(self)):
            return self.dot(other)
        else:
            raise TypeError(f"unsupported operand type(s) for *: '{type(self).__name__}' and '{type(other).__name__}'")

    def __matmul__(self, other):
        if isinstance(other, type(self)):
            return self.dot(other)
        else:
            raise TypeError(f"unsupported operand type(s) for @: '{type(self).__name__}' and '{type(other).__name__}'")

    def __add__(self, other):
        if isinstance(other, type(self)):
            res = Hammer.FourMomentum(self)
            res += other
            return res
        else:
            raise TypeError(f"unsupported operand type(s) for +: '{type(self).__name__}' and '{type(other).__name__}'")

    def __sub__(self, other):
        if isinstance(other, type(self)):
            res = Hammer.FourMomentum(self)
            res -= other
            return res
        else:
            raise TypeError(f"unsupported operand type(s) for -: '{type(self).__name__}' and '{type(other).__name__}'")

    def __truediv__(self, other):
        if isinstance(other, (float, int)):
            res = Hammer.FourMomentum(self)
            res /= other
            return res
        else:
            raise TypeError(f"unsupported operand type(s) for /: '{type(self).__name__}' and '{type(other).__name__}'")

    def __eq__(self, other):
        if self is other:
            return True
        return self.E() == other.E() and self.px() == other.px() and self.py() == other.py() and self.pz() == other.pz()

    if name == 'FourMomentum':
        klass.__len__ = lambda self: 4
        klass.__repr__ = lambda self: repr(str(self))
        klass.__str__ = __str__
        klass.__mul__ = __mul__
        klass.__rmul__ = __mul__
        klass.__add__ = __add__
        klass.__radd__ = __add__
        klass.__sub__ = __sub__
        klass.__matmul__ = __matmul__
        klass.__truediv__ = __truediv__
        klass.__eq__ = __eq__

add_pythonization(pythonize_four_momentum, 'Hammer')

def pythonize_particle(klass, name):
# turn pdgId and momentum into properties?
    def __str__(self):
        return f"Particle({self.p()}, {self.pdgId()})"

    if name == 'Particle':
        klass.__repr__ = lambda self: repr(str(self))
        klass.__str__ = __str__

add_pythonization(pythonize_particle, 'Hammer')
