import hammer

def toParticle(p):
    return hammer.Particle(hammer.FourMomentum(p.momentum().e(), p.momentum().px(), p.momentum().py(), p.momentum().pz()),
                            p.pdg_id())

def addDecayAndProducts(proc, p, parentId):
    daughters =[]
    for po in p.children():
        daughters.append(proc.addParticle(toParticle(po)));
        if po.status() == 2 and po.end_vertex() and po.pdg_id() != 111:
            addDecayAndProducts(proc, po, daughters[-1])
    proc.addVertex(parentId, daughters)

def parseGenEvent(evt, pdgCodes):
    results =[]
    for pi in evt.particles():
        if pi.status() == 2 and pi.pdg_id() in pdgCodes:
            p = hammer.Process()
            idxParent = p.addParticle(toParticle(pi))
            addDecayAndProducts(p, pi, idxParent)
            results.append(p)
    return results
