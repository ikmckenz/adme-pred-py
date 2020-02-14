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

    def druglikeness_lipinski(self, verbose=False):
        violations = []

        h_bond_donors = self.h_bond_donors()
        if h_bond_donors > 5:
            violations.append("H Bond Donors {}>5".format(h_bond_donors))

        h_bond_acceptors = self.h_bond_acceptors()
        if h_bond_acceptors > 10:
            violations.append("H Bond Acceptors {}>10".format(h_bond_acceptors))

        molecular_weight = self.molecular_weight()
        if molecular_weight > 500:
            violations.append("Molecular Weight {}>500".format(molecular_weight))

        logp = self.logp()
        if logp > 5:
            violations.append("LOGP {}>5".format(logp))

        if verbose:
            return violations
        else:
            return len(violations) < 1

    def druglikeness_ghose(self, verbose=False):
        violations = []

        logp = self.logp()
        if logp > 5.6 or logp < -0.4:
            violations.append("LOGP {}".format(logp))

        molecular_mass = self.molecular_weight()
        if molecular_mass < 160 or molecular_mass > 480:
            violations.append("Molecular Mass {}".format(molecular_mass))

        molar_refractivity = self.molar_refractivity()
        if molar_refractivity < 40 or molar_refractivity > 130:
            violations.append("MR {}".format(molar_refractivity))

        n_atoms = self.n_atoms()
        if n_atoms < 20 or n_atoms > 70:
            violations.append("N Atoms {}".format(n_atoms))

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

    def molar_refractivity(self):
        return Chem.Crippen.MolMR(self.mol)

    def molecular_weight(self):
        return Descriptors.ExactMolWt(self.mol)

    def n_atoms(self):
        return self.mol.GetNumAtoms()

    def logp(self):
        return Descriptors.MolLogP(self.mol)

    def tpsa(self):
        return Descriptors.TPSA(self.mol)


if __name__ == "__main__":
    chem = "ClC1=CC2=C(C=C1)N3C(C)=NN=C3CN=C2C4=CC=CC=C4"
    mol = ADME(chem)

    print(mol.druglikeness_ghose())
