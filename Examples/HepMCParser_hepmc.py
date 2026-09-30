import hammer
from hammer.hepmc import HepMCParser as _HepMCParser


class GenEvent:
    """Thin mutable container filled by ReaderAsciiHepMC2.read_event()."""
    def __init__(self):
        self._event = None


class ReaderAsciiHepMC2:
    def __init__(self, filename):
        self._parser = _HepMCParser(filename)
        self._failed = False

    def read_event(self, ge):
        event = self._parser.read_event()
        if event is None:
            self._failed = True
        else:
            ge._event = event

    def failed(self):
        return self._failed


def _add_decay(proc, pcode, parent_idx, evt):
    daughters = []
    for child_id in evt.particles[pcode].end_vertex(evt).children:
        child = evt.particles[child_id]
        idx = proc.addParticle(hammer.Particle(
            hammer.FourMomentum(child.p.e, child.p.px, child.p.py, child.p.pz),
            child.pdg
        ))
        daughters.append(idx)
        if child.status_code == 2 and child.end_vertex_code != 0 and child.pdg != 111:
            _add_decay(proc, child_id, idx, evt)
    proc.addVertex(parent_idx, daughters)


def parseGenEvent(ge, pdgCodes):
    if ge._event is None:
        return []
    evt = ge._event
    results = []
    for pcode, pi in evt.particles.items():
        if pi.status_code == 2 and pi.pdg in pdgCodes:
            p = hammer.Process()
            parent_idx = p.addParticle(hammer.Particle(
                hammer.FourMomentum(pi.p.e, pi.p.px, pi.p.py, pi.p.pz),
                pi.pdg
            ))
            _add_decay(p, pcode, parent_idx, evt)
            results.append(p)
    return results
