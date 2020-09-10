import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import FilterCatalog
from rdkit.Chem import rdqueries


class ADME(object):
    """
    ADME is the main class that contains the core ADME calculations. Functions
    include many druglikeness filters (Lipinski's rule of five, Egan, Ghose,
    BOILED-Egg, etc.), medicinal chemistry filters (PAINS, Brenk), and provides
    a consistent, simple wrapping API for common cheminformatics calculations.
    """

    # Class constants
    BOILED_EGG_HIA_ELLIPSE = Ellipse((71.051, 2.292), 142.081, 8.740, -1.031325)  # From the BOILED-Egg paper
    BOILED_EGG_BBB_ELLIPSE = Ellipse((38.117, 3.177), 82.061, 5.557, -0.171887)  # From the BOILED-Egg paper

    def __init__(self, mol):
        if isinstance(mol, str):
            self.mol = Chem.MolFromSmiles(mol)
        else:
            self.mol = mol

    def full_report(self):
        """Print a full report on a molecule with many implemented functions"""

        print("Pharmacokinetics:")
        print("GI absorption: {}".format("High" if self.boiled_egg_bbb() else "Low"))
        print("BBB permeant: {}".format("Yes" if self.boiled_egg_hia() else "No"))
        print(os.linesep)

        druglikeness = self.druglikeness_lipinski(verbose=True)
        print("Lipinski Rule of 5 Violations:")
        if isinstance(druglikeness, str):
            print(druglikeness)
        else:
            print(*druglikeness, sep="\n")
        print(os.linesep)

        print("Medicinal Chemistry:")
        print("PAINS filter: {}".format(self.pains()))
        print("Brenk filter: {}".format(self.brenk()))

    def druglikeness_egan(self, verbose=False):
        """
        Egan (2000) Prediction of Drug Absorption Using Multivariate Statistics

        The Egan paper creates the "Egan egg" a multivariable ellipse model.
        For simplicity, we don't implement the actual Egan egg here, but rather
        the simpler rule based approximation that the paper describes as the
        upper bounds for their model inputs. We use the 95% confidence interval
        levels here, rather than the 99% confidence intervals.

        In the future the full Egan egg model should be implemented here, but
        that is low priority given the problems with the Egan model.
        """
        violations = []

        tpsa = self.tpsa()
        if tpsa > 131.6:
            violations.append("PSA {}".format(tpsa))

        logp = self.logp()
        if logp > 5.88:
            violations.append("logP {}".format(logp))

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
        """
        Lipinski (2001) Experimental and computational approaches to estimate
        solubility and permeability in drug discovery and development settings

        https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five

        Lipinski's rule of 5 is one of the most important druglikenss filters,
        against which all others are judged. The rules of the filter are no
        more than 5 hydrogen bond donors, no more than 10 hydrogen bond
        acceptors, a molecular mass of less than 500 daltons, and a logP
        that does not exceed 5.
        """
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
            if not violations:
                return "No violations found"
            return violations
        else:
            return len(violations) < 1

    def druglikeness_muegge(self, verbose=False):
        """
        Muegge (2001) Simple Selection Criteria for Drug-like Chemical Matter

        In this paper they define a few specific "pharmacophore points" that
        are essentially proxies of hydrogen bonding ability. SwissADME uses
        a few simple rules based filters to approximate the pharmacophore point
        filter in this paper, and we implement the rules based filter here.
        """
        violations = []

        molecular_weight = self.molecular_weight()
        if molecular_weight > 600 or molecular_weight < 200:
            violations.append("MW {}".format(molecular_weight))

        logp = self.logp()
        if logp > 5 or logp < -2:
            violations.append("LOGP {}".format(logp))

        tpsa = self.tpsa()
        if tpsa > 150:
            violations.append("TPSA {}".format(tpsa))

        n_rings = self.n_rings()
        if n_rings > 7:
            violations.append("N Rings {}".format(n_rings))

        n_carbon = self.n_carbons()
        if n_carbon < 5:
            violations.append("N Carbon {}".format(n_carbon))

        n_heteroatoms = self.n_heteroatoms()
        if n_heteroatoms < 2:
            violations.append("N Heteroatoms".format(n_heteroatoms))

        n_rot_bonds = self.n_rot_bonds()
        if n_rot_bonds > 15:
            violations.append("N Rot Bonds {}".format(n_rot_bonds))

        h_bond_acc = self.h_bond_acceptors()
        if h_bond_acc > 10:
            violations.append("H Bond Acc {}".format(h_bond_acc))

        h_bond_don = self.h_bond_donors()
        if h_bond_don > 5:
            violations.append("H Bond Don {}".format(h_bond_don))

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

    def boiled_egg_bbb(self, logp=None, psa=None):
        """
        Daina (2016) A BOILED-Egg To Predict Gastrointestinal Absorption and
        Brain Penetration of Small Molecules

        This multivariate model uses log P and Polar Surface Area to determine
        druglikeness. This function implements their Blood Brain Barrier
        (BBB) model, which is the "yolk" of the BOILED-Egg.
        """

        if logp is None:
            logp = self.logp()

        if psa is None:
            psa = self.tpsa()

        return self.BOILED_EGG_BBB_ELLIPSE.contains_point((psa, logp))

    def boiled_egg_hia(self, logp=None, psa=None):
        """
        Daina (2016) A BOILED-Egg To Predict Gastrointestinal Absorption and
        Brain Penetration of Small Molecules

        This multivariate model uses log P and Polar Surface Area to determine
        druglikeness. This function implements their Human Intestinal Absorption
        (HIA) model, which is the "white" of the BOILED-Egg. 
        """

        if logp is None:
            logp = self.logp()

        if psa is None:
            psa = self.tpsa()

        return self.BOILED_EGG_HIA_ELLIPSE.contains_point((psa, logp))

    def boiled_egg_graphical(self):
        """
        Daina (2016) A BOILED-Egg To Predict Gastrointestinal Absorption and
        Brain Penetration of Small Molecules

        Takes the BOILED-Egg calculations implemented in boiled_egg_bbb and
        boiled_egg_hia and creates a matplotlib plot of the computations.
        :return:
        """
        fig, ax = plt.subplots()

        ax.patch.set_facecolor("lightgrey")

        white = self.BOILED_EGG_HIA_ELLIPSE
        yolk = self.BOILED_EGG_BBB_ELLIPSE

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

    def brenk(self):
        """
        Brenk (2008) Lessons Learnt from Assembling Screening Libraries for
        Drug Discovery for Neglected Diseases

        Brenk's Structural Alert filter finds fragments "putatively toxic,
        chemically reactive, metabolically unstable or to bear properties
        responsible for poor pharmacokinetics."

        Returns:
            Boolean of whether the molecule triggers the Brenk filter.
        """
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        catalog = FilterCatalog.FilterCatalog(params)
        return catalog.HasMatch(self.mol)

    def pains(self):
        """
        Baell and Holloway (2010) New Substructure Filters for Removal of Pan
        Assay Interference Compounds (PAINS) from Screening Libraries and for
        Their Exclusion in Bioassays

        This filter finds promiscuous compounds that are likely to show activity
        regardless of the target.

        Returns:
            Boolean of whether the molecule triggers the PAINS filter.
        """
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)
        return catalog.HasMatch(self.mol)

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

    def n_carbons(self):
        carbon = Chem.rdqueries.AtomNumEqualsQueryAtom(6)
        return len(self.mol.GetAtomsMatchingQuery(carbon))

    def n_heteroatoms(self):
        return Descriptors.rdMolDescriptors.CalcNumHeteroatoms(self.mol)

    def n_rings(self):
        return Descriptors.rdMolDescriptors.CalcNumRings(self.mol)

    def n_rot_bonds(self):
        return Chem.Lipinski.NumRotatableBonds(self.mol)

    def logp(self):
        return Descriptors.MolLogP(self.mol)

    def tpsa(self):
        return Descriptors.TPSA(self.mol)


if __name__ == "__main__":
    chem = "O=C(C)Oc1ccccc1C(=O)O"
    mol = ADME(chem)
    mol.full_report()