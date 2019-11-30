import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from rdkit import Chem
from rdkit.Chem import Descriptors


class ADME(object):

    def __init__(self, mol):
        if isinstance(mol, str):
            self.mol = Chem.MolFromSmiles(mol)
        else:
            self.mol = mol

    def boiled_egg(self):
        fig, ax = plt.subplots()

        ax.patch.set_facecolor("lightgrey")

        white = Ellipse((71.051, 2.292), 142.081, 8.740, -1.031325)
        yolk = Ellipse((38.117, 3.177), 82.061, 5.557, -0.171887)

        white.set_clip_box(ax.bbox)
        white.set_facecolor("white")
        ax.add_artist(white)

        yolk.set_clip_box(ax.bbox)
        yolk.set_facecolor("orange")
        ax.add_artist(yolk)

        ax.set_xlim(-10, 200)
        ax.set_ylim(-4, 8)

        ax.set_xlabel("TPSA")
        ax.set_ylabel("LOGP")

        logp = self.logp()
        tpsa = self.tpsa()

        ax.scatter(tpsa, logp, zorder=10)

        return fig

    def logp(self):
        return Descriptors.MolLogP(self.mol)

    def tpsa(self):
        return Descriptors.TPSA(self.mol)


if __name__ == "__main__":
    chem = "ClC1=CC2=C(C=C1)N3C(C)=NN=C3CN=C2C4=CC=CC=C4"
    mol = ADME(chem)

    fig = mol.boiled_egg()
    plt.show()
