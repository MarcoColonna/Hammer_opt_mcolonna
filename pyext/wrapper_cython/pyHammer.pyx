#**** This file is a part of the HAMMER library
#**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
#**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
#**** Please note the MCnet academic guidelines; see GUIDELINES for details

## \file  pyHammer.pyx
#  \brief Cython class to wrap Hammer library
#

# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=default

from typing import List, Set, Dict, Tuple, Optional, FrozenSet, AnyStr, Union, Any, NoReturn
from numbers import Number



Vec3 = Tuple[float, float, float]

from cppdefs cimport Hammer as cpp_Hammer
from cppdefs cimport Process as cpp_Process
from cppdefs cimport Particle as cpp_Particle
from cppdefs cimport FourMomentum as cpp_FourMomentum
from cppdefs cimport IOBuffer as cpp_IOBuffer
from cppdefs cimport IOBuffers as cpp_IOBuffers
from cppdefs cimport WTerm as cpp_WTerm
from cppdefs cimport PAction as cpp_PAction
from cppdefs cimport HistoInfo as cpp_HistoInfo
from cppdefs cimport RecordType as cpp_RecordType
from cppdefs cimport BinContents as cpp_BinContents
from cppdefs cimport UMap as cpp_UMap
from cppdefs cimport Log as cpp_Log
from cppdefs cimport version as cpp_version

from cpython cimport array
import array
from cpython cimport Py_buffer

from io import FileIO
from struct import pack, unpack, calcsize
from enum import Enum, IntEnum

from libc.stdint cimport uint8_t, uint16_t
from libcpp.string cimport string
from libcpp.map cimport map as cmap
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr
from libcpp.set cimport set as cset
from cython.operator import dereference as deref, preincrement as inc
from cpython cimport bool

from cymove cimport cymove as move
from cpython.version cimport PY_MAJOR_VERSION

import cmath

HAVE_NUMPY = True
try:
    import numpy as np
    import numpy.typing as npt
except:
    HAVE_NUMPY = False

def _is_string(value: Any) -> bool:
    return isinstance(value, (str, bytes, bytearray))

def _is_set(value: Any) -> bool:
    return isinstance(value, (set, frozenset))

def _is_event_uid_group(value: Any) -> bool:
    return (_is_set(value) and
            all([(_is_set(elem) and
                all([(type(elem2) is int and elem2>0)
                for elem2 in elem]))
            for elem in value]))

def _is_list_of_string(value: Any) -> bool:
    return (isinstance(value, list) and all([_is_string(elem) for elem in value]))

def _is_list_of_list_of_string(value: Any) -> bool:
    return (isinstance(value, list) and all([_is_list_of_string(elem) for elem in value]))

def _is_vec_float(value: Any) -> bool:
    return (isinstance(value, list) and all([type(elem) is float for elem in value]))

def _is_vec_complex(value: Any) -> bool:
    return (isinstance(value, list) and all([isinstance(elem, Number) for elem in value]))

def _is_vec_int(value: Any) -> bool:
    return (isinstance(value, list) and all([type(elem) is int for elem in value]))

def _is_dict_float(value: Any) -> bool:
    return (isinstance(value, dict) and all([(type(elem_v) is float and _is_string(elem_k)) for elem_k, elem_v in value.items()]))

def _is_dict_complex(value: Any) -> bool:
    return (isinstance(value, dict) and all([(isinstance(elem_v, Number) and _is_string(elem_k)) for elem_k, elem_v in value.items()]))

cdef _to_set_of_sets(cset[cset[size_t]] value):
    return frozenset([frozenset([int(elem2) for elem2 in elem]) for elem in value])

cdef class EventUIDGroup():

    cdef cset[cset[size_t]] c_data

    def __init__(self, ids):
        if ids is not None:
            self.__from_python(ids)

    @staticmethod
    cdef EventUIDGroup from_cpp(cset[cset[size_t]] ids):
        cdef EventUIDGroup wrapper = EventUIDGroup.__new__(EventUIDGroup)
        wrapper.c_data = ids
        return wrapper

    def __from_python(self, data):
        cdef cset[size_t] elem_out
        self.c_data.clear()
        try:
            for elem in data:
                elem_out = elem
                self.c_data.insert(elem_out)
                elem_out.clear()
        except:
            raise TypeError("Cannot convert python data into EventUIDGroup")

    cdef cset[cset[size_t]] to_cpp(self):
        return self.c_data
    
cdef class BinSizes():

    cdef vector[uint16_t] c_data

    def __init__(self, ids=None):
        if ids is not None:
            self.__from_python(ids)

    @staticmethod
    cdef BinSizes from_cpp(vector[uint16_t] ids):
        cdef BinSizes wrapper = BinSizes.__new__(BinSizes)
        wrapper.c_data = ids
        return wrapper

    def __from_python(self, data):
        try:
            self.c_data = data
        except:
            raise TypeError("Cannot convert python data into BinSizes")

    cdef vector[uint16_t] to_cpp(self):
        return self.c_data

cdef class BinEdges():

    cdef vector[vector[double]] c_data

    def __init__(self, ids=None):
        if ids is not None:
            self.__from_python(ids)

    @staticmethod
    cdef BinEdges from_cpp(vector[vector[double]] ids):
        cdef BinEdges wrapper = BinEdges.__new__(BinEdges)
        wrapper.c_data = ids
        return wrapper

    def __from_python(self, data):
        cdef vector[double] elem_out
        self.c_data.clear()
        try:
            for elem in data:
                elem_out = elem
                self.c_data.push_back(elem_out)
                elem_out.clear()
        except:
            raise TypeError("Cannot convert python data into BinEdges")

    cdef vector[vector[double]] to_cpp(self):
        return self.c_data

cdef class BinRanges():

    cdef vector[pair[double, double]] c_data

    def __init__(self, ids=None):
        if ids is not None:
            self.__from_python(ids)

    @staticmethod
    cdef BinRanges from_cpp(vector[pair[double, double]] ids):
        cdef BinRanges wrapper = BinRanges.__new__(BinRanges)
        wrapper.c_data = ids
        return wrapper

    def __from_python(self, data):
        cdef pair[double, double] elem_out
        self.c_data.clear()
        try:
            for elem in data:
                elem_out = elem
                self.c_data.push_back(elem_out)
        except:
            raise TypeError("Cannot convert python data into BinEdges")

    cdef vector[pair[double, double]] to_cpp(self):
        return self.c_data

cdef class IntList():

    cdef vector[int] c_data

    def __init__(self, ids):
        if ids is not None:
            self.__from_python(ids)

    @staticmethod
    cdef IntList from_cpp(vector[int] ids):
        cdef IntList wrapper = IntList.__new__(IntList)
        wrapper.c_data = ids
        return wrapper

    def __from_python(self, data):
        try:
            self.c_data = data
        except:
            raise TypeError("Cannot convert python data into IntList")

    cdef vector[int] to_cpp(self):
        return self.c_data

class WTerm(Enum):
    COMMON = 0
    NUMERATOR = 1
    DENOMINATOR = 2

cdef cpp_WTerm _to_cpp_wterm(arg):
    if arg == WTerm.COMMON:
        return cpp_WTerm.COMMON
    elif arg == WTerm.NUMERATOR:
        return cpp_WTerm.NUMERATOR
    elif arg == WTerm.DENOMINATOR:
        return cpp_WTerm.DENOMINATOR
    else:
        return cpp_WTerm.COMMON

class PAction(Enum):
    ALL = 0
    WEIGHTS = 1
    HISTOGRAMS = 2

cdef cpp_PAction _to_cpp_paction(arg):
    if arg == PAction.ALL:
        return cpp_PAction.ALL
    elif arg == PAction.WEIGHTS:
        return cpp_PAction.WEIGHTS
    elif arg == PAction.HISTOGRAMS:
        return cpp_PAction.HISTOGRAMS
    else:
        return cpp_PAction.ALL

class RecordType(IntEnum):
    UNDEFINED = ord('u')
    HEADER = ord('b')
    EVENT = ord('e')
    HISTOGRAM = ord('h')
    RATE = ord('r')
    HISTOGRAM_DEFINITION = ord('d')

cdef class HistoInfo:

    cdef cpp_HistoInfo c_histo_info

    cdef from_cpp(self, cpp_HistoInfo info):
        self.c_histo_info = info

    # Attribute access
    @property
    def name(self) -> AnyStr:
        return self.c_histo_info.name
    @name.setter
    def name(self, name: AnyStr) -> NoReturn:
        self.c_histo_info.name = name

    @property
    def scheme(self) -> AnyStr:
        return self.c_histo_info.scheme
    @scheme.setter
    def scheme(self, scheme: AnyStr) -> NoReturn:
        self.c_histo_info.scheme = scheme

    @property
    def event_group_id(self) -> EventUIDGroup:
        return set([frozenset(elem) for elem in self.c_histo_info.eventGroupId])
    @event_group_id.setter
    def event_group_id(self, event_group_id: EventUIDGroup) -> NoReturn:
        self.c_histo_info.eventGroupId = event_group_id

cdef class BinContents:

    cdef cpp_BinContents c_bins

    @staticmethod
    cdef from_cpp(cpp_BinContents value):
        result = BinContents()
        result.c_bins = value
        return result

    @property
    def sum_wi(self) -> float:
        return self.c_bins.sumWi
    @sum_wi.setter
    def sum_wi(self, sum_wi: float) -> NoReturn:
        self.c_bins.sumWi = sum_wi

    @property
    def sum_wi2(self) -> float:
        return self.c_bins.sumWi2
    @sum_wi2.setter
    def sum_wi2(self, sum_wi2: float) -> NoReturn:
        self.c_bins.sumWi2 = sum_wi2

    @property
    def n(self) -> int:
        return self.c_bins.n
    @n.setter
    def n(self, n: int) -> NoReturn:
        self.c_bins.n = n

cdef class IOBuffers:

    cdef cpp_IOBuffers* p_io_buffers
    cdef cpp_IOBuffers.iterator curr_it

    def __cinit__(self):
        self.p_io_buffers = new cpp_IOBuffers()
        self.curr_it = self.p_io_buffers.begin()

    def __dealloc__(self):
        del self.p_io_buffers

    cdef from_cpp(self, cpp_IOBuffers buf):
        del self.p_io_buffers
        self.p_io_buffers = new cpp_IOBuffers(move(buf))
        self.curr_it = self.p_io_buffers.begin()

    def __iter__(self):
        self.curr_it = self.p_io_buffers.begin()
        return self

    def __next__(self):
        if self.curr_it == self.p_io_buffers.end():
            raise StopIteration
        result = IOBuffer()
        result.from_cpp(deref(self.curr_it))
        inc(self.curr_it)
        return result

    def save(self, file_handler : FileIO) -> NoReturn:
        cdef cpp_IOBuffers.iterator it = self.p_io_buffers.begin()
        while it != self.p_io_buffers.end():
            result = IOBuffer()
            result.from_cpp(deref(it))
            result.save(file_handler)
            inc(it)

cdef class IOBuffer:

    cdef cpp_IOBuffer c_buffer
    cdef array.array buf
    cdef int view_count
    cdef Py_ssize_t strides

    def __cinit__(self):
        self.c_buffer.kind = cpp_RecordType.UNDEFINED
        self.c_buffer.length = 0
        self.c_buffer.start = NULL
        self.view_count = 0

    def __len__(self):
        return self.c_buffer.length

    def __getbuffer__(self, Py_buffer *view, int flags):
        cdef Py_ssize_t itemsize = 1
        view.buf = <unsigned char *> self.c_buffer.start
        view.format = 'B'
        view.internal = NULL
        view.itemsize = 1
        view.len = self.c_buffer.length
        view.ndim = 1
        view.obj = self
        view.readonly = 0
        view.shape = <Py_ssize_t*> &self.c_buffer.length
        view.strides = &itemsize
        view.suboffsets = NULL
        self.view_count += 1

    def __releasebuffer__(self, Py_buffer *buffer):
        self.view_count -= 1

    cdef from_cpp(self, cpp_IOBuffer buf):
        if self.view_count > 0:
            raise ValueError("can't edit buffer while it's accessed via buffer protocol")
        self.c_buffer = buf
        self.view_count = 0

    cdef cpp_IOBuffer* to_cpp(self):
        return &self.c_buffer

    cpdef init_from_size(self, int data_size = 0):
        if self.view_count > 0:
            raise ValueError("can't edit buffer while it's accessed via buffer protocol")
        self.c_buffer.init(data_size)
        self.view_count = 0

    @property
    def kind(self) -> RecordType:
        return self.c_buffer.kind

    @kind.setter
    def kind(self, kind: RecordType) -> NoReturn:
        self.c_buffer.kind = <cpp_RecordType>(kind)

    def save(self, file_handler) -> NoReturn:
        file_handler.write(pack(b'b', self.c_buffer.kind))
        file_handler.write(pack(b'<L', self.c_buffer.length))
        file_handler.write(self)

    def load(self, file_handler) -> bool:
        tmp = file_handler.read(1)
        if len(tmp) == 0:
            return False
        typebuf = unpack(b'b', tmp)
        form = b'<L'
        sizebuf = unpack(form, file_handler.read(calcsize(form)))[0]
        self.init_from_size(sizebuf)
        self.c_buffer.kind = <cpp_RecordType>(typebuf[0])
        file_handler.readinto(self)
        return True


cdef class FourMomentum:

    cdef cpp_FourMomentum c_mom

    def __cinit__(self, e: float = 0., px: float = 0., py: float = 0., pz: float = 0.):
        self.c_mom = cpp_FourMomentum(e, px, py, pz)

    @staticmethod
    def fromPtEtaPhiM(pt: float, eta: float, phi: float, m: float) -> FourMomentum:
        result = FourMomentum()
        result.c_mom = cpp_FourMomentum.fromPtEtaPhiM(pt, eta, phi, m)
        return result

    @staticmethod
    def fromEtaPhiME(eta: float, phi: float, m: float, e: float) -> FourMomentum:
        result = FourMomentum()
        result.c_mom = cpp_FourMomentum.fromEtaPhiME(eta, phi, m, e)
        return result

    @staticmethod
    def fromPM(px: float, py: float, pz: float, m: float) -> FourMomentum:
        result = FourMomentum()
        result.c_mom = cpp_FourMomentum.fromPM(px, py, pz, m)
        return result

    @property
    def px(self) -> float:
        return self.c_mom.px()
    @px.setter
    def px(self, value: float) -> NoReturn:
        self.c_mom.setPx(value)

    @property
    def py(self) -> float:
        return self.c_mom.py()
    @py.setter
    def py(self, value: float) -> NoReturn:
        self.c_mom.setPy(value)

    @property
    def pz(self) -> float:
        return self.c_mom.pz()
    @pz.setter
    def pz(self, value: float) -> NoReturn:
        self.c_mom.setPz(value)

    @property
    def e(self) -> float:
        return self.c_mom.E()
    @e.setter
    def e(self, value: float) -> NoReturn:
        self.c_mom.setE(value)

    def mass(self) -> float:
        return self.c_mom.mass()

    def mass2(self) -> float:
        return self.c_mom.mass2()

    def p(self) -> float:
        return self.c_mom.p()

    def p2(self) -> float:
        return self.c_mom.p2()

    def pt(self) -> float:
        return self.c_mom.pt()

    def rapidity(self) -> float:
        return self.c_mom.rapidity()

    def phi(self) -> float:
        return self.c_mom.phi()

    def eta(self) -> float:
        return self.c_mom.eta()

    def theta(self) -> float:
        return self.c_mom.theta()

    def p_vec(self) -> Vec3:
        res = self.c_mom.pVec()
        return (res[0], res[1], res[2])

    def gamma(self) -> float:
        return self.c_mom.gamma()

    def beta(self) -> float:
        return self.c_mom.beta()

    def boost_vector(self) -> Vec3:
        res = self.c_mom.boostVector()
        return (res[0], res[1], res[2])

    def dot(self, other: FourMomentum) -> float:
        return self.c_mom.dot((<FourMomentum>other).c_mom)

    def __len__(self) -> int:
        return 4

    def __str__(self) -> str:
        return f"FourMomentum({self.c_mom.E()}, {self.c_mom.px()}, {self.c_mom.py()}, {self.c_mom.pz()})"

    def __repr__(self) -> str:
        return repr(str(self))

    def __add__(self, other):
        if isinstance(other, FourMomentum):
            return FourMomentum(self.c_mom.E() + (<FourMomentum>other).c_mom.E(),
                                self.c_mom.px() + (<FourMomentum>other).c_mom.px(),
                                self.c_mom.py() + (<FourMomentum>other).c_mom.py(),
                                self.c_mom.pz() + (<FourMomentum>other).c_mom.pz())
        raise TypeError(f"unsupported operand type(s) for +: 'FourMomentum' and '{type(other).__name__}'")

    def __sub__(self, other):
        if isinstance(other, FourMomentum):
            return FourMomentum(self.c_mom.E() - (<FourMomentum>other).c_mom.E(),
                                self.c_mom.px() - (<FourMomentum>other).c_mom.px(),
                                self.c_mom.py() - (<FourMomentum>other).c_mom.py(),
                                self.c_mom.pz() - (<FourMomentum>other).c_mom.pz())
        raise TypeError(f"unsupported operand type(s) for -: 'FourMomentum' and '{type(other).__name__}'")

    def __mul__(self, other):
        if isinstance(other, (float, int)):
            return FourMomentum(self.c_mom.E() * other, self.c_mom.px() * other,
                                self.c_mom.py() * other, self.c_mom.pz() * other)
        elif isinstance(other, FourMomentum):
            return self.c_mom.dot((<FourMomentum>other).c_mom)
        raise TypeError(f"unsupported operand type(s) for *: 'FourMomentum' and '{type(other).__name__}'")

    def __rmul__(self, other):
        if isinstance(other, (float, int)):
            return FourMomentum(self.c_mom.E() * other, self.c_mom.px() * other,
                                self.c_mom.py() * other, self.c_mom.pz() * other)
        raise TypeError(f"unsupported operand type(s) for *: '{type(other).__name__}' and 'FourMomentum'")

    def __truediv__(self, other):
        if isinstance(other, (float, int)):
            return FourMomentum(self.c_mom.E() / other, self.c_mom.px() / other,
                                self.c_mom.py() / other, self.c_mom.pz() / other)
        raise TypeError(f"unsupported operand type(s) for /: 'FourMomentum' and '{type(other).__name__}'")

    def __matmul__(self, other):
        if isinstance(other, FourMomentum):
            return self.c_mom.dot((<FourMomentum>other).c_mom)
        raise TypeError(f"unsupported operand type(s) for @: 'FourMomentum' and '{type(other).__name__}'")

    def __eq__(self, other) -> bool:
        if not isinstance(other, FourMomentum):
            return NotImplemented
        return (self.c_mom.E() == (<FourMomentum>other).c_mom.E() and
                self.c_mom.px() == (<FourMomentum>other).c_mom.px() and
                self.c_mom.py() == (<FourMomentum>other).c_mom.py() and
                self.c_mom.pz() == (<FourMomentum>other).c_mom.pz())

cdef class Particle:

    cdef cpp_Particle c_part

    def __cinit__(self, p: FourMomentum = FourMomentum(), pdg: int = 0):
        self.c_part = cpp_Particle(p.c_mom, pdg)

    @staticmethod
    cdef Particle from_cpp(cpp_Particle p):
        part = Particle()
        part.c_part = p
        return part

    @property
    def momentum(self) -> FourMomentum:
        result = FourMomentum()
        result.c_mom = self.c_part.momentum()
        return result

    @momentum.setter
    def momentum(self, val: FourMomentum) -> NoReturn:
        self.c_part.setMomentum(val.c_mom)

    @property
    def p(self) -> FourMomentum:
        result = FourMomentum()
        result.c_mom = self.c_part.momentum()
        return result

    @p.setter
    def p(self, val: FourMomentum) -> NoReturn:
        self.c_part.setMomentum(val.c_mom)

    @property
    def pdg_id(self) -> int:
        return self.c_part.pdgId()

    @pdg_id.setter
    def pdg_id(self, val: int) -> NoReturn:
        self.c_part.setPdgId(val)

    def __str__(self) -> str:
        return f"Particle({self.momentum}, {self.c_part.pdgId()})"

    def __repr__(self) -> str:
        return repr(str(self))


cdef class Process:

    cdef cpp_Process c_proc

    def add_particle(self, value: Particle) -> NoReturn:
        return self.c_proc.addParticle(value.c_part)

    def add_vertex(self, parent: int, daughters: List[int]) -> NoReturn:
        self.c_proc.addVertex(parent, daughters)

    def remove_vertex(self, vertex_id: int, prune: bool = False) -> NoReturn:
        self.c_proc.removeVertex(vertex_id, prune)

    def get_daughters_ids(self, parent: int = 0) -> List[int]:
        return self.c_proc.getDaughtersIds(parent)

    def get_parent_id(self, daughter: int) -> int:
        return self.c_proc.getParentId(daughter)

    def get_siblings_ids(self, particle: int = 0) -> List[int]:
        return self.c_proc.getSiblingsIds(particle)

    def get_particle(self, id: int) -> Particle:
        result = Particle()
        result.c_part = self.c_proc.getParticle(id)
        return result

    def get_daughters(self, parent: int = 0, sorted: bool = False) -> List[Particle]:
        result = self.c_proc.getDaughters(parent, sorted)
        return [Particle.from_cpp(elem) for elem in result]

    def get_siblings(self, particle: int = 0, sorted: bool = False) -> List[Particle]:
        result = self.c_proc.getSiblings(particle, sorted)
        return [Particle.from_cpp(elem) for elem in result]

    def get_parent(self, daughter: int) -> Particle:
        result = Particle()
        result.c_part = self.c_proc.getParent(daughter)
        return result

    def is_parent(self, particle: int) -> bool:
        return self.c_proc.isParent(particle)

    def get_first_vertex(self) -> int:
        return self.c_proc.getFirstVertex()

    def get_id(self) -> int:
        return self.c_proc.getId()

    def full_id(self) -> Set[int]:
        return self.c_proc.fullId()

    def get_vertex_id(self, vertex) -> int:
        return self.c_proc.getVertexId(vertex)

    def num_particles(self, without_photons: bool = False) -> int:
        return self.c_proc.numParticles(without_photons)

    def get_particles_by_vertex(self, parent_or_vertex, daughters=None) -> Tuple[Particle, List[Particle]]:
        if daughters is not None:
            result = self.c_proc.getParticlesByVertex(parent_or_vertex, daughters)
        else:
            result = self.c_proc.getParticlesByVertex(parent_or_vertex)
        particles = [Particle.from_cpp(item) for item in result.second]
        return (Particle.from_cpp(result.first), particles)



## \brief main Hammer python class
#
#  Provides the functionalities of `Hammer::Hammer` visible to a python program
#
#  \ingroup PyExt

cdef class Hammer:

    ## pointer to Hammer::Hammer class
    cdef cpp_Hammer *wrapped

    ## base constructor
    def __cinit__ (self):
        self.wrapped = new cpp_Hammer()
        if self.wrapped is NULL:
            raise MemoryError()

    ## destructor
    def __dealloc__ (self):
        if self.wrapped is not NULL:
            del self.wrapped

    def init_run(self) -> NoReturn:
        self.wrapped.initRun()

    def init_event(self, weight: float = 1.0) -> NoReturn:
        self.wrapped.initEvent(weight)

    def add_process(self, proc: Process) -> int:
        return self.wrapped.addProcess(proc.c_proc)

    def remove_process(self, proc_id: int) -> NoReturn:
        self.wrapped.removeProcess(proc_id)

    def set_event_histogram_bin(self, name: AnyStr, bins: List[int]) -> NoReturn:
        self.wrapped.setEventHistogramBin(name, bins)

    def fill_event_histogram(self, name: AnyStr, values: List[float]) -> NoReturn:
        self.wrapped.fillEventHistogram(name, values)

    def set_event_base_weight(self, weight: float) -> NoReturn:
        self.wrapped.setEventBaseWeight(weight)

    def process_event(self, what: PAction = PAction.ALL) -> NoReturn:
        self.wrapped.processEvent(_to_cpp_paction(what))

    def load_event_weights(self, buf: IOBuffer, merge: bool = False) -> bool:
        return self.wrapped.loadEventWeights(deref(buf.to_cpp()), merge)

    def save_event_weights(self) -> IOBuffer:
        result = IOBuffer()
        result.from_cpp(self.wrapped.saveEventWeights())
        return result

    def write_event_weights(self, file_handler: FileIO) -> NoReturn:
        file_handler.write(self.save_event_weights())

    def load_run_header(self, buf: IOBuffer, merge: bool = False) -> bool:
        return self.wrapped.loadRunHeader(deref(buf.to_cpp()), merge)

    def save_run_header(self) -> IOBuffer:
        result = IOBuffer()
        result.from_cpp(self.wrapped.saveRunHeader())
        return result

    def write_run_header(self, file_handler: FileIO) -> NoReturn:
        file_handler.write(self.save_run_header())

    def load_histogram_definition(self, buf: IOBuffer, merge: bool = False) -> AnyStr:
        return self.wrapped.loadHistogramDefinition(deref(buf.to_cpp()), merge)

    def load_histogram(self, buf: IOBuffer, merge: bool = False) -> HistoInfo:
        result = HistoInfo()
        result.from_cpp(self.wrapped.loadHistogram(deref(buf.to_cpp()), merge))
        return result

    def save_histogram(self, *args) -> IOBuffers:
        cdef string s_name
        cdef string s_scheme
        cdef string s_spec
        cdef cset[cset[size_t]] ids
        result = IOBuffers()
        if not args:
            raise TypeError("save_histogram requires at least one argument")
        if isinstance(args[0], HistoInfo):
            result.from_cpp(self.wrapped.saveHistogram((<HistoInfo>args[0]).c_histo_info))
        elif len(args) == 1:
            s_name = args[0]
            result.from_cpp(self.wrapped.saveHistogram(s_name))
        elif len(args) == 2 and isinstance(args[1], EventUIDGroup):
            s_name = args[0]
            ids = (<EventUIDGroup>args[1]).to_cpp()
            result.from_cpp(self.wrapped.saveHistogram(s_name, ids))
        elif len(args) == 2:
            s_name = args[0]
            s_scheme = args[1]
            result.from_cpp(self.wrapped.saveHistogram(s_name, s_scheme))
        elif len(args) == 3 and isinstance(args[2], EventUIDGroup):
            s_name = args[0]
            s_scheme = args[1]
            ids = (<EventUIDGroup>args[2]).to_cpp()
            result.from_cpp(self.wrapped.saveHistogram(s_name, s_scheme, ids))
        elif len(args) == 3:
            s_name = args[0]
            s_scheme = args[1]
            s_spec = args[2]
            result.from_cpp(self.wrapped.saveHistogram(s_name, s_scheme, s_spec))
        elif len(args) == 4:
            s_name = args[0]
            s_scheme = args[1]
            s_spec = args[2]
            ids = (<EventUIDGroup>args[3]).to_cpp()
            result.from_cpp(self.wrapped.saveHistogram(s_name, s_scheme, s_spec, ids))
        return result

#    def _save_histo_namedimpl(self, **kwargs) -> IOBuffers:
#        result = IOBuffers()
#        cdef string name
#        cdef string scheme
#        cdef cset[cset[size_t]] ids
#        if 'info' in kwargs:
#            result.from_cpp(self.wrapped.saveHistogram(((HistoInfo)(kwargs['info'])).c_histo_info))
#        elif 'name' in kwargs:
#            name = kwargs['name']
#            if 'scheme' in  kwargs:
#                scheme = kwargs['scheme']
#                if 'event_ids' in kwargs:
#                    ids = kwargs['event_ids']
#                    result.from_cpp(self.wrapped.saveHistogram(name, scheme, ids))
#                else:
#                    result.from_cpp(self.wrapped.saveHistogram(name, scheme))
#            else:
#                if 'event_ids' in kwargs:
#                    ids = kwargs['event_ids']
#                    result.from_cpp(self.wrapped.saveHistogram(name, ids))
#                else:
#                    result.from_cpp(self.wrapped.saveHistogram(name))
#        else:
#            raise NotImplementedError("Invalid arguments")
#        return result

    def write_histogram(self, file_handler: FileIO, *args, **kwargs) -> NoReturn:
        for elem in self.save_histogram(*args, **kwargs):
            file_handler.write(elem)

    def load_rates(self, buf: IOBuffer, merge: bool = False) -> bool:
        return self.wrapped.loadRates(deref(buf.to_cpp()), merge)

    def save_rates(self) -> IOBuffer:
        result = IOBuffer()
        result.from_cpp(self.wrapped.saveRates())
        return result

    def write_rates(self, file_handler: FileIO) -> NoReturn:
        file_handler.write(self.save_rates())

    def read_cards(self, file_decays: AnyStr, file_options: AnyStr) -> NoReturn:
        self.wrapped.readCards(file_decays, file_options)

    def save_option_card(self, file_options: AnyStr, use_default: bool = True) -> NoReturn:
        self.wrapped.saveOptionCard(file_options, use_default)

    def save_header_card(self, file_decays: AnyStr) -> NoReturn:
        self.wrapped.saveHeaderCard(file_decays)

    def save_references(self, file_refs: AnyStr) -> NoReturn:
        self.wrapped.saveReferences(file_refs)

    def set_options(self, options: AnyStr) -> NoReturn:
        self.wrapped.setOptions(options)

    def set_header(self, options: AnyStr) -> NoReturn:
        self.wrapped.setHeader(options)

    def add_total_sum_of_weights(self, compress: bool = False, with_errors: bool = False) -> NoReturn:
        self.wrapped.addTotalSumOfWeights(compress, with_errors)

    def add_histogram(self, name, bins, has_under_over_flow=False, ranges=None) -> NoReturn:
        cdef string s_name = name
        if isinstance(bins, BinSizes):
            if ranges is None:
                ranges = BinRanges()
            self.wrapped.addHistogram(s_name, (<BinSizes>bins).to_cpp(), has_under_over_flow, (<BinRanges>ranges).to_cpp())
        elif isinstance(bins, BinEdges):
            self.wrapped.addHistogram(s_name, (<BinEdges>bins).to_cpp(), has_under_over_flow)
        else:
            raise TypeError(f"Expected BinSizes or BinEdges, got {type(bins)}")

#    def __add_histogram_namedimpl(self, **kwargs):
#        if 'name' not in kwargs:
#            raise NotImplementedError("Invalid arguments")
#        cdef string name = kwargs['name']
#        has_under_over_flow = kwargs.get('has_under_over_flow', True)
#        if 'bin_sizes' in kwargs:
#            ranges = kwargs.get('ranges',[])
#            self.wrapped.addHistogram(name, kwargs['bin_sizes'].to_cpp(), has_under_over_flow, ranges)
#        elif 'bin_edges' in kwargs:
#            self.wrapped.addHistogram(name, kwargs['bin_edges'].to_cpp(), has_under_over_flow)
#        else:
#            raise NotImplementedError("Invalid arguments")

    def collapse_processes_in_histogram(self, name: AnyStr) -> NoReturn:
        self.wrapped.collapseProcessesInHistogram(name)

    def keep_errors_in_histogram(self, name: AnyStr, value: bool = True) -> NoReturn:
        self.wrapped.keepErrorsInHistogram(name, value)

    def reconcile_specializations(self) -> NoReturn:
        self.wrapped.reconcileSpecializations()

    def save_event_general_weights(self, save: bool = True) -> NoReturn:
        self.wrapped.saveEventGeneralWeights(save)

    def reconcile_general_weight_setting(self) -> NoReturn:
        self.wrapped.reconcileGeneralWeightSetting()

    def reconcile_general_histogram_setting(self, histogram_name: AnyStr) -> NoReturn:
        self.wrapped.reconcileGeneralHistogramSetting(histogram_name)

    def create_wc_specialization(self, name: AnyStr, wc_space: AnyStr, coordinates: List[AnyStr] = []) -> NoReturn:
        self.wrapped.createWCSpecialization(name, wc_space, coordinates)

    def set_wc_specialization_origin(self, name: AnyStr, wc_space: AnyStr, values_or_settings: Union[List[complex], Dict[AnyStr, complex]]) -> NoReturn:
        cdef vector[double complex] v
        cdef cmap[string, double complex] d
        cdef string s_name = name
        cdef string s_space = wc_space
        if _is_vec_complex(values_or_settings):
            v = values_or_settings
            self.wrapped.setWCSpecializationOrigin(s_name, s_space, v)
        elif _is_dict_complex(values_or_settings):
            d = values_or_settings
            self.wrapped.setWCSpecializationOrigin(s_name, s_space, d)
        else:
            print("Invalid Type")

    def set_wc_specialization_coord(self, name: AnyStr, wc_space: AnyStr, coord: AnyStr, values_or_settings: Union[List[complex], Dict[AnyStr, complex]]) -> NoReturn:
        cdef vector[double complex] v
        cdef cmap[string, double complex] d
        cdef string s_name = name
        cdef string s_space = wc_space
        cdef string s_coord = coord
        if _is_vec_complex(values_or_settings):
            v = values_or_settings
            self.wrapped.setWCSpecializationCoord(s_name, s_space, s_coord, v)
        elif _is_dict_complex(values_or_settings):
            d = values_or_settings
            self.wrapped.setWCSpecializationCoord(s_name, s_space, s_coord, d)
        else:
            print("Invalid Type")

    def set_wc_specialization_basis(self, name: AnyStr, wc_space: AnyStr, sub_space: List[Dict[AnyStr, complex]]) -> NoReturn:
        self.wrapped.setWCSpecializationBasis(name, wc_space, sub_space)

    def remove_wc_specialization(self, name: AnyStr, remove_existing_data: bool = False) -> NoReturn:
        self.wrapped.removeWCSpecialization(name, remove_existing_data)

    def remove_all_wc_specializations(self, remove_existing_data: bool = False) -> NoReturn:
        self.wrapped.removeAllWCSpecializations(remove_existing_data)

    def available_wc_specialization_ids(self) -> Set[AnyStr]:
        return self.wrapped.availableWCSpecializationIds()

    def apply_wc_specialization_in_weights(self, name: AnyStr) -> NoReturn:
        self.wrapped.applyWCSpecializationInWeights(name)

    def remove_wc_specialization_in_weights(self, name: AnyStr) -> NoReturn:
        self.wrapped.removeWCSpecializationInWeights(name)

    def applied_wc_specializations_in_weights(self) -> Set[AnyStr]:
        return self.wrapped.appliedWCSpecializationsInWeights()

    def apply_wc_specialization_in_histogram(self, histogram_name: AnyStr, spec_name: AnyStr) -> NoReturn:
        self.wrapped.applyWCSpecializationInHistogram(histogram_name, spec_name)

    def apply_wc_specialization_in_all_histograms(self, spec_name: AnyStr) -> NoReturn:
        self.wrapped.applyWCSpecializationInAllHistograms(spec_name)

    def remove_wc_specialization_in_histogram(self, histogram_name: AnyStr, spec_name: AnyStr) -> NoReturn:
        self.wrapped.removeWCSpecializationInHistogram(histogram_name, spec_name)

    def applied_wc_specializations_in_histogram(self, histogram_name: AnyStr) -> Set[AnyStr]:
        return self.wrapped.appliedWCSpecializationsInHistogram(histogram_name)

#    def specialize_wc_in_weights(self, process: AnyStr, values_or_settings: Union[List[complex], Dict[AnyStr, complex]]) -> NoReturn:
#        cdef string s1 = process
#        cdef vector[double complex] v
#        cdef cmap[string, double complex] d
#        if _is_vec_complex(values_or_settings):
#            v = values_or_settings
#            self.wrapped.specializeWCInWeights(s1, v)
#        elif _is_dict_complex(values_or_settings):
#            d = values_or_settings
#            self.wrapped.specializeWCInWeights(s1, d)
#        else:
#            print("Invalid Type")

 #   def specialize_wc_in_histogram(self, name: AnyStr, wc_space: AnyStr, values_or_settings: Union[List[complex], Dict[AnyStr, complex]]) -> NoReturn:
 #       cdef vector[double complex] v
 #       cdef cmap[string, double complex] d
 #       cdef string s_name = name
 #       cdef string s_space = wc_space
 #       if _is_vec_complex(values_or_settings):
 #           v = values_or_settings
 #           self.wrapped.specializeWCInHistogram(s_name, s_space, v)
 #       elif _is_dict_complex(values_or_settings):
 #           d = values_or_settings
 #           self.wrapped.specializeWCInHistogram(s_name, s_space, d)
 #       else:
 #           print("Invalid Type")

    def specialize_ff_in_histogram(self, name: AnyStr, process: AnyStr, group: AnyStr, values_or_settings: Union[List[float], Dict[AnyStr, float]]) -> NoReturn:
        cdef vector[double] v
        cdef cmap[string, double] d
        cdef string s_name = name
        cdef string s_process = process
        cdef string s_group = group
        if _is_vec_float(values_or_settings):
            v = values_or_settings
            self.wrapped.specializeFFInHistogram(s_name, s_process, s_group, v)
        elif _is_dict_float(values_or_settings):
            d = values_or_settings
            self.wrapped.specializeFFInHistogram(s_name, s_process, s_group, d)
        else:
            print("Invalid Type")

#    def reset_specialization_in_histogram(self, name: AnyStr) -> NoReturn:
#        self.wrapped.resetSpecializationInHistogram(name)

    def remove_ff_specialization_in_histogram(self, name: AnyStr) -> NoReturn:
        self.wrapped.removeFFSpecializationInHistogram(name)

    def create_projected_histogram(self, old_name: AnyStr, new_name: AnyStr, collapsed_index_positions: Set[int]) -> NoReturn:
        self.wrapped.createProjectedHistogram(old_name, new_name, collapsed_index_positions)

    def remove_histogram(self, name: AnyStr) -> NoReturn:
        self.wrapped.removeHistogram(name)

    def add_ff_scheme(self, scheme_name: AnyStr, schemes: Dict[AnyStr, AnyStr]) -> NoReturn:
        self.wrapped.addFFScheme(scheme_name, schemes)

    def set_ff_input_scheme(self, schemes: Dict[AnyStr, AnyStr]) -> NoReturn:
        self.wrapped.setFFInputScheme(schemes)

    def remove_ff_scheme(self, scheme_name: AnyStr) -> NoReturn:
        self.wrapped.removeFFScheme(scheme_name)

    def get_ff_scheme_names(self) -> List[str]:
        res = self.wrapped.getFFSchemeNames()
        return [str(elem) for elem in res]

    def show_available_ff_params(self, prefix: AnyStr = '') -> NoReturn:
        cdef string s
        if prefix:
            s = prefix
            self.wrapped.showAvailableFFParams(s)
        else:
            self.wrapped.showAvailableFFParams()

    def include_decay(self, name_or_names: Union[List[AnyStr], AnyStr]) -> NoReturn:
        cdef string s
        cdef vector[string] v
        if _is_list_of_string(name_or_names):
            v = name_or_names
            self.wrapped.includeDecay(v)
        elif _is_string(name_or_names):
            s = name_or_names
            self.wrapped.includeDecay(s)
        else:
            print("Invalid Type")

    def forbid_decay(self, name_or_names: Union[List[AnyStr], AnyStr]) -> NoReturn:
        cdef string s
        cdef vector[string] v
        if _is_list_of_string(name_or_names):
            v = name_or_names
            self.wrapped.forbidDecay(v)
        elif _is_string(name_or_names):
            s = name_or_names
            self.wrapped.forbidDecay(s)
        else:
            print("Invalid Type")

    def add_pure_ps_vertices(self, vertices: Set[AnyStr], what: WTerm = WTerm.NUMERATOR) -> NoReturn:
        self.wrapped.addPurePSVertices(vertices, _to_cpp_wterm(what))

    def clear_pure_ps_vertices(self, what: WTerm) -> NoReturn:
        self.wrapped.clearPurePSVertices(_to_cpp_wterm(what))

    def set_units(self, name: AnyStr = "GeV") -> NoReturn:
        self.wrapped.setUnits(name)

    def rename_ff_eigenvectors(self, process: AnyStr, group: AnyStr, names: List[AnyStr]) -> NoReturn:
        cdef vector[string] v
        v = names
        self.wrapped.renameFFEigenvectors(process, group, v)

    def set_wilson_coefficients(self, wc_space: AnyStr, values_or_settings: Union[List[complex], Dict[AnyStr, complex]], what = WTerm.NUMERATOR) -> NoReturn:
        cdef vector[double complex] v
        cdef cmap[string, double complex] d
        cdef string s_space = wc_space
        if _is_vec_complex(values_or_settings):
            v = values_or_settings
            self.wrapped.setWilsonCoefficients(s_space, v, _to_cpp_wterm(what))
        elif _is_dict_complex(values_or_settings):
            d = values_or_settings
            self.wrapped.setWilsonCoefficients(s_space, d, _to_cpp_wterm(what))
        else:
            print("Invalid Type")

    def set_wilson_coefficients_local(self, wc_space: AnyStr, values_or_settings: Union[List[complex], Dict[AnyStr, complex]]) -> NoReturn:
        cdef vector[double complex] v
        cdef cmap[string, double complex] d
        cdef string s_space = wc_space
        if _is_vec_complex(values_or_settings):
            v = values_or_settings
            self.wrapped.setWilsonCoefficientsLocal(s_space, v)
        elif _is_dict_complex(values_or_settings):
            d = values_or_settings
            self.wrapped.setWilsonCoefficientsLocal(s_space, d)
        else:
            print("Invalid Type")

    def reset_wilson_coefficients(self, wc_space: AnyStr, what: WTerm = WTerm.NUMERATOR) -> NoReturn:
        self.wrapped.resetWilsonCoefficients(wc_space, _to_cpp_wterm(what))

    def retrieve_wilson_coefficients(self, wc_space: AnyStr, what: WTerm = WTerm.NUMERATOR) -> Dict[AnyStr, complex]:
        cdef string s_space = wc_space
        return self.wrapped.retrieveWilsonCoefficients(s_space, _to_cpp_wterm(what))

    def set_ff_eigenvectors(self, process: AnyStr, group: AnyStr, values_or_settings: Union[List[float], Dict[AnyStr, float]]) -> NoReturn:
        cdef vector[double] v
        cdef cmap[string, double] d
        cdef string s_process = process
        cdef string s_group = group
        if _is_vec_float(values_or_settings):
            v = values_or_settings
            self.wrapped.setFFEigenvectors(s_process, s_group, v)
        elif _is_dict_float(values_or_settings):
            d = values_or_settings
            self.wrapped.setFFEigenvectors(s_process, s_group, d)
        else:
            print("Invalid Type")

    def set_ff_eigenvectors_local(self, process: AnyStr, group: AnyStr, values_or_settings: Union[List[float], Dict[AnyStr, float]]) -> NoReturn:
        cdef vector[double] v
        cdef cmap[string, double] d
        cdef string s_process = process
        cdef string s_group = group
        if _is_vec_float(values_or_settings):
            v = values_or_settings
            self.wrapped.setFFEigenvectorsLocal(s_process, s_group, v)
        elif _is_dict_float(values_or_settings):
            d = values_or_settings
            self.wrapped.setFFEigenvectorsLocal(s_process, s_group, d)
        else:
            print("Invalid Type")

    def reset_ff_eigenvectors(self, process: AnyStr, group: AnyStr) -> NoReturn:
        self.wrapped.resetFFEigenvectors(process, group)

    def retrieve_ff_eigenvectors(self, process: AnyStr, group: AnyStr) -> Dict[AnyStr, float]:
        cdef string s_process = process
        cdef string s_group = group
        return self.wrapped.retrieveFFEigenvectors(s_process, s_group)

    def get_weight(self, scheme, processes_or_specialization=None, specialization='') -> float:
        cdef string s_scheme = scheme
        cdef string s_spec
        cdef vector[size_t] v1
        cdef vector[vector[string]] v2
        if processes_or_specialization is None or _is_string(processes_or_specialization):
            s_spec = processes_or_specialization if processes_or_specialization is not None else ''
            return self.wrapped.getWeight(s_scheme, s_spec)
        else:
            s_spec = specialization
            if _is_vec_int(processes_or_specialization):
                v1 = processes_or_specialization
                return self.wrapped.getWeight(s_scheme, v1, s_spec)
            elif _is_list_of_list_of_string(processes_or_specialization):
                v2 = processes_or_specialization
                return self.wrapped.getWeight(s_scheme, v2, s_spec)
            else:
                print("Invalid Type")
                return 0.

    def get_weights(self, scheme: AnyStr, specialization: AnyStr = '') -> Dict[int, float]:
        return self.wrapped.getWeights(scheme, specialization)


    def get_rate(self, *args) -> float:
        cdef Py_ssize_t i_id
        cdef int i_parent
        cdef string s_vertex
        cdef string s_scheme
        cdef string s_spec
        if not args or len(args) < 2:
            raise TypeError("get_rate requires at least two arguments")
        if isinstance(args[0], str):
            s_vertex = args[0]
            s_scheme = args[1]
            s_spec = args[2] if len(args) > 2 else ''
            return self.wrapped.getRate(s_vertex, s_scheme, s_spec)
        elif len(args) >= 3 and isinstance(args[1], IntList):
            i_parent = args[0]
            s_scheme = args[2]
            s_spec = args[3] if len(args) > 3 else ''
            return self.wrapped.getRate(i_parent, (<IntList>args[1]).to_cpp(), s_scheme, s_spec)
        else:
            i_id = args[0]
            s_scheme = args[1]
            s_spec = args[2] if len(args) > 2 else ''
            return self.wrapped.getRate(i_id, s_scheme, s_spec)

#    def __get_rate_namedimpl(self, **kwargs):
#        if 'scheme' not in kwargs:
#            raise NotImplementedError("Invalid arguments")
#        cdef string scheme = kwargs['scheme']
#        cdef Py_ssize_t id
#        cdef string vertex
#        if 'id' in kwargs:
#            id = kwargs['id']
#            return self.wrapped.getRate(id,scheme)
#        elif 'vertex' in kwargs:
#            vertex = kwargs['vertex']
#            return self.wrapped.getRate(vertex,scheme)
#        elif ('parent' in kwargs) and ('daughters' in kwargs):
#            return self.wrapped.getRate(kwargs['parent'], kwargs['daughters'], scheme)
#        else:
#            raise NotImplementedError("Invalid arguments")

    def get_denominator_rate(self, *args) -> float:
        cdef Py_ssize_t i1
        cdef string s_vertex
        if not args:
            raise TypeError("get_denominator_rate requires at least one argument")
        if isinstance(args[0], str):
            s_vertex = args[0]
            return self.wrapped.getDenominatorRate(s_vertex)
        elif len(args) >= 2 and isinstance(args[1], IntList):
            return self.wrapped.getDenominatorRate(args[0], (<IntList>args[1]).to_cpp())
        else:
            i1 = args[0]
            return self.wrapped.getDenominatorRate(i1)

#    def _get_denominator_rate_namedimpl(self, **kwargs):
#        cdef Py_ssize_t id
#        cdef string vertex
#        if 'id' in kwargs:
#            id = kwargs['id']
#            return self.wrapped.getDenominatorRate(id)
#        elif 'vertex' in kwargs:
#            vertex = kwargs['vertex']
#            return self.wrapped.getDenominatorRate(vertex)
#        elif ('parent' in kwargs) and ('daughters' in kwargs):
#            return self.wrapped.getDenominatorRate(kwargs['parent'], kwargs['daughters'])
#        else:
#            raise NotImplementedError("Invalid arguments")

    def get_histogram(self, name: AnyStr, scheme: AnyStr, specialization: AnyStr = '') -> List[BinContents]:
        result = self.wrapped.getHistogram(name, scheme, specialization)
        return [BinContents.from_cpp(elem) for elem in result]

    def get_histograms(self, name: AnyStr, scheme: AnyStr, specialization: AnyStr = '') -> Dict[FrozenSet[FrozenSet[int]], List[BinContents]]:
        cdef cpp_UMap[cset[cset[size_t]], vector[cpp_BinContents]] res
        cdef cpp_UMap[cset[cset[size_t]], vector[cpp_BinContents]].iterator it
        res = self.wrapped.getHistograms(name, scheme, specialization)
        results = dict()
        it = res.begin()
        while (it != res.end()):
            tmp_v = [BinContents.from_cpp(elem) for elem in deref(it).second]
            tmp_k = _to_set_of_sets(deref(it).first)
            results[tmp_k] = tmp_v
            inc(it)
        return results

    if HAVE_NUMPY:

        def get_numpy_histogram(self, name: AnyStr, scheme: AnyStr, specialization: AnyStr = '') -> npt.NDArray[BinContents]:
            shape = self.wrapped.getHistogramShape(name)
            res = self.wrapped.getHistogram(name, scheme, specialization)
            results = np.array([BinContents.from_cpp(o) for o in res], dtype=BinContents)
            results.reshape(list(shape))
            return results

        def get_numpy_histograms(self, name: AnyStr, scheme: AnyStr, specialization: AnyStr = '') -> Dict[FrozenSet[FrozenSet[int]], npt.NDArray[BinContents]]:
            cdef cpp_UMap[cset[cset[size_t]], vector[cpp_BinContents]] res
            cdef cpp_UMap[cset[cset[size_t]], vector[cpp_BinContents]].iterator it
            shape = self.wrapped.getHistogramShape(name)
            res = self.wrapped.getHistograms(name, scheme, specialization)
            results = dict()
            it = res.begin()
            while (it != res.end()):
                tmp_v = np.array([BinContents.from_cpp(o) for o in deref(it).second], dtype=BinContents)
                tmp_k = _to_set_of_sets(deref(it).first)
                results[tmp_k] = tmp_v
                inc(it)
            return results

    def get_histogram_event_ids(self, name: AnyStr, scheme: AnyStr, specialization: AnyStr = '') -> EventUIDGroup:
        return self.wrapped.getHistogramEventIds(name, scheme, specialization)

    def get_histogram_bin_edges(self, name: AnyStr) -> List[List[float]]:
        return self.wrapped.getHistogramBinEdges(name)

    def get_histogram_shape(self, name: AnyStr) -> List[int]:
        return self.wrapped.getHistogramShape(name)

    def histogram_has_under_over_flows(self, name: AnyStr) -> bool:
        return self.wrapped.histogramHasUnderOverFlows(name)

def version() -> str:
    return str(cpp_version())

class LogLevel(Enum):
    TRACE = 0
    DEBUG = 10
    INFO = 20
    WARN = 30
    WARNING = 30
    ERROR = 40
    CRITICAL = 50
    ALWAYS = 50

def set_log_level(name: AnyStr, level: LogLevel) -> NoReturn:
    cpp_Log.setLevel(name, int(level))

def set_log_levels(level_dict: Dict[AnyStr, LogLevel]) -> NoReturn:
    cpp_Log.setLevels(level_dict)

def set_log_show_timestamp(value: bool) -> NoReturn:
    cpp_Log.setShowTimestamp(value)

def set_log_show_level(value: bool) -> NoReturn:
    cpp_Log.setShowLevel(value)

def set_log_show_logger_name(value: bool) -> NoReturn:
    cpp_Log.setShowLoggerName(value)

def set_log_use_colors(value: bool) -> NoReturn:
    cpp_Log.setUseColors(value)

def set_log_warning_max_count(name: AnyStr, max_count: int) -> NoReturn:
    cpp_Log.setWarningMaxCount(name, max_count)

def reset_log_warning_counters() -> NoReturn:
    cpp_Log.resetWarningCounters()
