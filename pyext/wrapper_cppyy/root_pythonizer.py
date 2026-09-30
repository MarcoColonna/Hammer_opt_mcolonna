#**** This file is a part of the HAMMER library
#**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
#**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
#**** Please note the MCnet academic guidelines; see GUIDELINES for details


from io import BufferedIOBase
from struct import calcsize, pack, unpack

from cppyy import ll
from ROOT import Hammer, pythonization  # type: ignore[import-not-found]


@pythonization('FourMomentum', ns='Hammer')
def pythonize_four_momentum(klass):
    def _str__(self):
        return f"FourMomentum({self.E()}, {self.px()}, {self.py()}, {self.pz()})"
    
    def _mul__(self, other):
        if isinstance(other, (float, int)):
            res = Hammer.FourMomentum(self)
            res *= other
            return res
        elif isinstance(other, type(self)):
            return self.dot(other)
        else:
            raise TypeError(f"unsupported operand type(s) for *: '{type(self).__name__}' and '{type(other).__name__}'")

    def _matmul__(self, other):
        if isinstance(other, type(self)):
            res = Hammer.FourMomentum(self)
            return res.dot(other)
        else:
            raise TypeError(f"unsupported operand type(s) for @: '{type(self).__name__}' and '{type(other).__name__}'")

    def _add__(self, other):
        if isinstance(other, type(self)):
            res = Hammer.FourMomentum(self)
            res += other
            return res
        else:
            raise TypeError(f"unsupported operand type(s) for +: '{type(self).__name__}' and '{type(other).__name__}'")

    def _sub__(self, other):
        if isinstance(other, type(self)):
            res = Hammer.FourMomentum(self)
            res -= other
            return res
        else:
            raise TypeError(f"unsupported operand type(s) for -: '{type(self).__name__}' and '{type(other).__name__}'")

    def _div__(self, other):
        if isinstance(other, (float, int)):
            res = Hammer.FourMomentum(self)
            res /= other
            return res
        else:
            raise TypeError(f"unsupported operand type(s) for /: '{type(self).__name__}' and '{type(other).__name__}'")

    def _eq__(self, other):
        if self is other:
            return True
        return self.E() == other.E() and self.px() == other.px() and self.py() == other.py() and self.pz() == other.pz()

    klass.__len__ = lambda self: 4
    klass.__str__ = _str__
    klass.__mul__ = _mul__
    klass.__rmul__ = _mul__
    klass.__add__ = _add__
    klass.__radd__ = _add__
    klass.__sub__ = _sub__
    klass.__matmul__ = _matmul__
    klass.__truediv__ = _div__
    klass.__eq__ = _eq__
    klass.__repr__ = lambda self: repr(_str__(self))

@pythonization('Particle', ns='Hammer')
def pythonize_particle(klass):
# turn pdgId and momentum into properties?
    def _str__(self):
        return f"Particle({self.p()}, {self.pdgId()})"
    
    klass.__str__ = _str__
    klass.__repr__ = lambda self: repr(_str__(self))

@pythonization('IOBuffer', ns='Hammer')
def pythonize_iobuffer(klass):
#turn kind, length, start into properties?

    def _len__(self):
        return self.length

    def _full_load(self, iostream):
        if isinstance(iostream, BufferedIOBase):
            return self._load_from_file(iostream)
        else:
            self._load_from_buffer(iostream)
            return self.length > 0

    def _full_save(self, iostream):
        if isinstance(iostream, BufferedIOBase):
            self._save_to_file(iostream)
        else:
            self._save_to_buffer(iostream)

    def save_to_file(self: Hammer.IOBuffer, file_handler: BufferedIOBase):
        file_handler.write(pack(b'c', bytes(self.kind[0],'utf-8')))
        file_handler.write(pack(b'<L', self.length))
        data = ll.cast['uint8_t*'](self.start)
        data.reshape((self.length,))
        file_handler.write(data)

    def load_from_file(self, file_handler: BufferedIOBase) -> bool:
        tmp = file_handler.read(1)
        if len(tmp) == 0:
            return False
        typebuf = unpack(b'c',tmp)
        form = b'<L'
        sizebuf = unpack(form, file_handler.read(calcsize(form)))[0]
        self.init(sizebuf)
        self.kind = Hammer.RecordType(typebuf[0])
        data = ll.cast['uint8_t*'](self.start)
        data.reshape((self.length,))
        file_handler.readinto(data)
        return True

    klass.__len__ = _len__
    klass._save_to_buffer = klass.save
    klass._load_from_buffer = klass.load
    klass._save_to_file = save_to_file
    klass._load_from_file = load_from_file
    klass.load = _full_load
    klass.save = _full_save


@pythonization('IOBuffers', ns='Hammer')
def pythonize_iobuffers(klass):

    def _len__(self):
        return self.size()

    def _full_save(self, iostream):
        if isinstance(iostream, BufferedIOBase):
            self._save_to_file(iostream)
        else:
            self._save_to_buffer(iostream)

    def save_to_file(self: Hammer.IOBuffers, file_handler: BufferedIOBase):
        for item in self:
            item.save(file_handler)

    klass.__len__ = _len__
    klass._save_to_buffer = klass.save
    klass._save_to_file = save_to_file
    klass.save = _full_save
