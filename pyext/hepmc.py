#**** This file is a part of the HAMMER library
#**** Copyright (C) 2016 - 2026 The HAMMER Collaboration
#**** HAMMER is licensed under version 3 of the GPL; see COPYING for details
#**** Please note the MCnet academic guidelines; see GUIDELINES for details

"""
@package hammer
@file hepmc.py
@brief Define a simple HepMC text parser
"""

from math import sqrt


class FourPosition:
    def __init__(self, x=0.0, y=0.0, z=0.0, t=0.0):
        self.t = t
        self.x = x
        self.y = y
        self.z = z

    def distance(self, other=None):
        if other is None:
            return sqrt(self.x**2. + self.y**2. + self.z**2.)
        return sqrt((other.x - self.x)**2. + (other.y - self.y)**2. + (other.z - self.z)**2.)

    def to_vector(self):
        return [self.t, self.x, self.y, self.z]

    def __sub__(self, other):
        return FourPosition(self.x - other.x, self.y - other.y, self.z - other.z, self.t - other.t)

    def __add__(self, other):
        return FourPosition(self.x + other.x, self.y + other.y, self.z + other.z, self.t + other.t)


class FourMomentum:
    def __init__(self, px=0.0, py=0.0, pz=0.0, e=0.0):
        self.e = e
        self.px = px
        self.py = py
        self.pz = pz

    def mass2(self):
        return self.e**2 - self.px**2 - self.py**2 - self.pz**2

    def mass(self):
        return sqrt(self.e**2 - self.px**2 - self.py**2 - self.pz**2)

    def to_vector(self):
        return [self.e, self.px, self.py, self.pz]

    def __sub__(self, other):
        return FourMomentum(self.px - other.px, self.py - other.py, self.pz - other.pz, self.e - other.e)

    def __add__(self, other):
        return FourMomentum(self.px + other.px, self.py + other.py, self.pz + other.pz, self.e + other.e)

    def p(self):
        return sqrt(self.px**2 + self.py**2 + self.pz**2)


class Particle:
    def __init__(self, pdg=0, p=None, status_code=0, start_vertex=0, end_vertex=0):
        self.pdg = pdg
        self.status_code = status_code
        self.p = p if p is not None else FourMomentum()
        self.start_vertex_code = start_vertex
        self.end_vertex_code = end_vertex

    def start_vertex(self, event=None):
        return event.vertices[self.start_vertex_code] if event is not None else None

    def end_vertex(self, event=None):
        if event is None:
            return None
        if self.end_vertex_code in event.vertices:
            return event.vertices[self.end_vertex_code]
        return Vertex()


class Vertex:
    def __init__(self, pos=None, parents=None, children=None):
        self.pos = pos if pos is not None else FourPosition()
        self.parents = set(p for p in parents if p > 0) if parents is not None else set()
        self.children = set(c for c in children if c > 0) if children is not None else set()

    def add_parents(self, parents=None):
        if parents:
            self.parents.update(p for p in parents if p > 0)

    def add_children(self, children=None):
        if children:
            self.children.update(c for c in children if c > 0)

    def in_particles(self, event=None):
        if event:
            return [event.particles[bar_code] for bar_code in self.parents]
        return []

    def out_particles(self, event=None):
        if event:
            return [event.particles[bar_code] for bar_code in self.children]
        return []

    def distance(self, other=None):
        return self.pos.distance(other.pos) if other is not None else 0.


class Event:
    def __init__(self, number=0, process_id=0, weights=None):
        self.number = number
        self.weights = weights if weights is not None else []
        self.process_id = process_id
        self.units: list = []
        self.cross_section: list = []
        self.particles: dict = {}
        self.vertices: dict = {}


class HepMCParser:
    def __init__(self, filename):
        self.version = None
        self._file = open(filename)
        self._line: str = ""
        status = True
        while status:
            status = self._read_hepmc_line()
            if status and self._line.startswith("HepMC::Version"):
                self.version = self._line.split()[1]
                break
        while status and not self._line.startswith("HepMC::IO_GenEvent-START_EVENT_LISTING"):
            self._read_hepmc_line()

    def _read_hepmc_line(self) -> bool:
        line = self._file.readline()
        if not line:
            return False
        self._line = line.strip(" \n")
        return True

    def read_event(self):
        if not self._line or self._line.startswith("HepMC::IO_GenEvent-END_EVENT_LISTING"):
            return None
        while self._line.startswith("HepMC"):
            self._read_hepmc_line()
        if not self._line.startswith("E "):
            return None
        number, process_id, num_vertices, weights = self._parse_event_line()
        event = Event(number, process_id, weights)
        while self._read_hepmc_line() and not self._line.startswith("V "):
            if self._line.startswith("U "):
                event.units = self._parse_units()
            elif self._line.startswith("C "):
                event.cross_section = self._parse_cross_section()
        in_particles: dict = {}
        for _ in range(num_vertices):
            v_code, vertex, num_in_particles, num_out_particles = self._parse_vertex()
            in_list = []
            out_list = []
            for _ in range(num_in_particles):
                if self._read_hepmc_line():
                    p_code, particle, out_v_code = self._parse_particle()
                    event.particles[p_code] = particle
                    in_list.append(p_code)
            vertex.add_parents(in_list)
            for _ in range(num_out_particles):
                if self._read_hepmc_line():
                    p_code, particle, out_v_code = self._parse_particle(v_code)
                    event.particles[p_code] = particle
                    out_list.append(p_code)
                    if out_v_code != 0:
                        in_particles.setdefault(out_v_code, []).append(p_code)
            vertex.add_children(out_list)
            event.vertices[v_code] = vertex
            self._read_hepmc_line()
        for key, value in in_particles.items():
            event.vertices[key].add_parents(value)
        return event

    def _parse_cross_section(self):
        tokens = self._line.split()
        return [float(x) for x in tokens[1:3]]

    def _parse_event_line(self):
        tokens = self._line.split()
        num_seeds = int(tokens[11])
        return int(tokens[1]), int(tokens[6]), int(tokens[8]), [float(x) for x in tokens[13 + num_seeds:]]

    def _parse_units(self):
        tokens = self._line.split()
        return tokens[1:3]

    def _parse_vertex(self):
        tokens = self._line.split()
        code = -int(tokens[1])
        pos = [float(x) for x in tokens[3:7]]
        v = Vertex(FourPosition(*pos))
        return code, v, int(tokens[7]), int(tokens[8])

    def _parse_particle(self, start_vertex_code=0):
        tokens = self._line.split()
        mom = [float(x) for x in tokens[3:7]]
        end_v = -int(tokens[11])
        p = Particle(int(tokens[2]), FourMomentum(*mom), int(tokens[8]), start_vertex_code, end_v)
        return int(tokens[1]), p, end_v
