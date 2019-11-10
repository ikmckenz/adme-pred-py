from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import seaborn as sns
sns.set()


class ADME(object):

    def __init__(self, mol):
        if isinstance(mol, str):
            self.mol = Chem.MolFromSmiles(mol)
        else:
            self.mol = mol

    def logp(self):
        return Descriptors.MolLogP(self.mol)

    def tpsa(self):
        return Descriptors.TPSA(self.mol)


chem = "CN(CCC1=CNC2=C1C=CC=C2)C"
mol = ADME(chem)

fig, ax = plt.subplots()

white = Ellipse((71.051, 2.292), 142.081, 8.740, -1.031325)
yolk = Ellipse((38.117, 3.177), 82.061, 5.557, -0.171887)

ax.add_artist(white)
white.set_clip_box(ax.bbox)
white.set_facecolor("white")

ax.add_artist(yolk)
yolk.set_clip_box(ax.bbox)
yolk.set_facecolor("orange")

ax.set_xlim(-10, 200)
ax.set_ylim(-4, 8)

logp = mol.logp()
tpsa = mol.tpsa()

ax.plot(tpsa, logp)

plt.show()
