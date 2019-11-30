# adme-pred-py

### Getting started
To get up and running, first create the environment:
```bash
conda create --name adme-pred-py
conda activate adme-pred-py
conda install -c conda-forge -n adme-pred-py rdkit
conda install -c conda-forge -n adme-pred-py matplotlib
```

We want to replicate most of the functionality of [SwissADME](https://www.nature.com/articles/srep42717). 
We want the core code to be abstract enough that it can produce a extensive report on a single compound like SwissADME, or be used as a filter in a screen. (Ex. `for compound in compounds if lipinski_druglikeness then ...`)
We will probably also want to implement some functionality from [Open Source Bayesian Models](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4478615/).
