#include <vector>
#include <set>

#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAsciiHepMC2.h"

#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include "Hammer/Process.hh"


inline Hammer::Particle toParticle(const HepMC3::ConstGenParticlePtr p) {
    return Hammer::Particle{{p->momentum().e(), p->momentum().px(), p->momentum().py(), p->momentum().pz()},
                            p->pdg_id()};
}

static inline void addDecayAndProducts(Hammer::Process& proc, const HepMC3::ConstGenParticlePtr p, size_t parent) {
    std::vector<size_t> daughters;
    // auto endv = p->end_vertex();
    for (const auto& po : p->children()) {
        daughters.push_back(proc.addParticle(toParticle(po)));
        if (po->status() == 2 && po->end_vertex() && po->pdg_id() != 111) {
            addDecayAndProducts(proc, po, daughters.back());
        }
    }
    proc.addVertex(parent, daughters);
}

static inline std::vector<Hammer::Process> parseGenEvent(const HepMC3::GenEvent& evt, std::set<int> pdgCodes) {
    std::vector<Hammer::Process> results;
    for (const auto& pi : evt.particles()) {
        if (pi->status() == 2 && pdgCodes.find(pi->pdg_id()) != pdgCodes.end()) {
            Hammer::Process p;
            addDecayAndProducts(p, pi, p.addParticle(toParticle(pi)));
            results.push_back(p);
        }
    }
    return results;
}
