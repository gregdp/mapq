
import _multiscale
import chimera
import mapq
import mapq.mapq
from CGLutil.AdaptiveTree import AdaptiveTree


print "Calc Q scores:"

from VolumeViewer import Volume
vols = chimera.openModels.list(modelTypes = [Volume])
if len(vols) == 0 :
    print " - no volumes loaded"
    exit(0)
dmap = vols[0]
print " - volume: %s" % dmap.name


from chimera import Molecule
mols = chimera.openModels.list(modelTypes = [Molecule])
if len(mols) == 0 :
    print " - no molecules loaded"
    exit(0)
mol = mols[0]
print " - model: %s" % mol.name
mapq.mapq.SetBBAts ( mol )

ats = [at for at in mol.atoms if not at.element.name == "H"]
points = _multiscale.get_atom_coordinates ( ats, transformed = False )
print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)
#allAtTree = None


#q, qcc = CalcQp ( mol, None, dmap, allAtTree=allAtTree )
q, qcc = mapq.mapq.CalcQp ( mol, None, dmap, allAtTree=allAtTree )

if q >= -1 :
    print ""
    print "Avg. Q score: %.3f" % q
    print ""
