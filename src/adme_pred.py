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

    def lipinski_druglikeness(self, verbose=False):
        violations = []

        h_bond_donors = self.h_bond_donors()
        if h_bond_donors > 5:
            violations.append("H Bond Donors {}>5".format(h_bond_donors))

        h_bond_acceptors = self.h_bond_acceptors()
        if h_bond_acceptors > 10:
            violations.append("H Bond Acceptors {}>10".format(h_bond_acceptors))

        molecular_mass = self.molecular_mass()
        if molecular_mass > 500:
            violations.append("Molecular Mass {}>500".format(molecular_mass))

        logp = self.logp()
        if logp > 5:
            violations.append("LOGP {}>5".format(logp))

        if verbose:
            return violations
        else:
            return len(violations) < 1

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

    def h_bond_donors(self):
        return Chem.Lipinski.NumHDonors(self.mol)

    def h_bond_acceptors(self):
        return Chem.Lipinski.NumHAcceptors(self.mol)

    def molecular_mass(self):
        return Descriptors.ExactMolWt(self.mol)

    def logp(self):
        return Descriptors.MolLogP(self.mol)

    def tpsa(self):
        return Descriptors.TPSA(self.mol)


if __name__ == "__main__":
    chem = "ClC1=CC2=C(C=C1)N3C(C)=NN=C3CN=C2C4=CC=CC=C4"
    mol = ADME(chem)

    print(mol.lipinski_druglikeness())
