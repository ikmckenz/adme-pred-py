# adme-pred-py

This library supports computational drug discovery by implementing several druglikenss filters, medicinal chemistry filters, and provides an easy to use wrapping API for common cheminformatics calculations.

We want to replicate most of the functionality of [SwissADME](https://www.nature.com/articles/srep42717). 
We want the core code to be abstract enough that it can produce a extensive report on a single compound like SwissADME, or be used as a filter in a screen. (Ex. `for compound in compounds if lipinski_druglikeness then ...`)
We will probably also want to implement some functionality from [Open Source Bayesian Models](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4478615/).

### Getting started
We use Anaconda as the base Python, so first install it from here: https://www.anaconda.com/products/individual.

To get up and running, first create the environment:
```bash
conda create --name adme-pred-py
conda activate adme-pred-py
conda install -c conda-forge -n adme-pred-py rdkit
conda install -c conda-forge -n adme-pred-py matplotlib
```

### Examples

Asprin is druglike because it does not violate the Lipinski Rule of 5:
```python
chem = "O=C(C)Oc1ccccc1C(=O)O"
mol = ADME(chem)
print(mol.druglikeness_lipinski(verbose=True))
```
```bash
No violations found
```

Paromomycin is not a small molecule drug, and so druglikeness filters will screen it out:
```python
chem = "O=S(=O)(O)O.O([C@H]3[C@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@H]1O[C@@H](CN)[C@@H](O)[C@H](O)[C@H]1N)[C@H]2O)[C@@H](O)[C@H](N)C[C@@H]3N)[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4N)CO"
mol = ADME(chem)
print(*mol.druglikeness_lipinski(verbose=True), sep="\n")
```
```bash
H Bond Donors 15>5
H Bond Acceptors 21>10
Molecular Weight 713.263680664>500
```

Thalidomide triggers the Brenk medicinal chemistry filter:
```python
chem = "O=C(N1C2CCC(NC2=O)=O)C3=CC=CC=C3C1=O"
mol = ADME(chem)
print(mol.brenk())
```
```bash
True
```

