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

    def druglikeness_egan(self, verbose=False):
        """
        Egan (2000) Prediction of Drug Absorption Using Multivariate Statistics


        """
        violations = []

        """Egan uses an ellipse? Not a box. So what is SwissADME using? Paper isn't explicit, so I sent an email to ask."""

        if verbose:
            return violations
        else:
            return len(violations) < 1

    def druglikeness_ghose(self, verbose=False):
        """
        Ghose (1999) A Knowledge-Based Approach in Designing Combinatorial or
        Medicinal Chemistry Libraries for Drug Discovery.
        1. A Qualitative and Quantitative Characterization of Known Drug
        Databases

        In this paper (and this function) they define the qualifying range
        (the chance of missing good compounds is less than 20%) as a certain
        range of Log P, molecular weight, molar refractivity, and the total
        number of atoms. They also define a preferred range as the  interval
        having 50% of drugs. The tighter preferred range is implemented in
        druglikeness_ghose_pref.
        """
        violations = []

        logp = self.logp()
        if logp > 5.6 or logp < -0.4:
            violations.append("LOGP {}".format(logp))

        molecular_weight = self.molecular_weight()
        if molecular_weight < 160 or molecular_weight > 480:
            violations.append("Molecular Mass {}".format(molecular_weight))

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

    def druglikeness_ghose_pref(self, verbose=False):
        """
        Ghose (1999) A Knowledge-Based Approach in Designing Combinatorial or
        Medicinal Chemistry Libraries for Drug Discovery.
        1. A Qualitative and Quantitative Characterization of Known Drug
        Databases

        In this paper they define the qualifying range (the chance of missing
        good compounds is less than 20%) as a certain range of Log P, molecular
        weight, molar refractivity, and the total number of atoms. They define
        a preferred range (implemented in this function) as the  interval
        having 50% of drugs. The looser qualifying range is implemented in
        druglikeness_ghose.
        """
        violations = []

        logp = self.logp()
        if logp > 4.1 or logp < 1.3:
            violations.append("LOGP {}".format(logp))

        molecular_weight = self.molecular_weight()
        if molecular_weight < 230 or molecular_weight > 390:
            violations.append("Molecular Mass {}".format(molecular_weight))

        molar_refractivity = self.molar_refractivity()
        if molar_refractivity < 70 or molar_refractivity > 110:
            violations.append("MR {}".format(molar_refractivity))

        n_atoms = self.n_atoms()
        if n_atoms < 30 or n_atoms > 55:
            violations.append("N Atoms {}".format(n_atoms))

        if verbose:
            return violations
        else:
            return len(violations) < 1

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


    def druglikeness_veber(self, verbose=False):
        """
        Veber (2002) Molecular Properties That Influence the Oral
        Bioavailability of Drug Candidates

        This study on oral bioavailability in rats shows that molecular weight
        is only a good predictor of bioavailability insofar as it's correlated
        with rotatable bonds and polar surface area. They actually propose two
        different druglikeness rules that achieve a similar effect:
        i. polar surface area <= 140 angstroms squared AND # rot. bonds <= 10
        OR
        ii. sum of H-bond donors and acceptors <= 12 AND # rot. bonds <= 10

        SwissADME uses the polar surface area metric
        """
        violations = []

        tpsa = self.tpsa()
        if tpsa > 140:
            violations.append("TPSA {}".format(tpsa))

        n_rot_bonds = self.n_rot_bonds()
        if n_rot_bonds > 10:
            violations.append("N Rotatable Bonds {}".format(n_rot_bonds))

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

    def n_rot_bonds(self):
        return Chem.Lipinski.NumRotatableBonds(self.mol)

    def logp(self):
        return Descriptors.MolLogP(self.mol)

    def tpsa(self):
        return Descriptors.TPSA(self.mol)


if __name__ == "__main__":
    chem = "ClC1=CC2=C(C=C1)N3C(C)=NN=C3CN=C2C4=CC=CC=C4"
    mol = ADME(chem)

    print(mol.druglikeness_veber())
