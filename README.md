# Molecular Visualiser

### Package Requirements:
-----------
 - Python 3.7+ or PyCharm 2019+
 - [Anaconda](https://www.anaconda.com/products/individual)
 - [RDKit](https://www.rdkit.org/docs/index.html) 
 - [Open-Source PyMOL](https://pymolwiki.org/index.php/Windows_Install)
 - [ChEMBL webresource client](https://chembl.gitbook.io/chembl-interface-documentation/web-services)
 - [Biopython](https://biopython.org/)
 
 These modules were installed under a conda environment: RDKit, Open-Source PyMOL, and Biopython. 
 'pip install' was used to install ChEMBL webresource client.

### How to Run:
-----------
To run the software workflow under interactive mode, open the Command Prompt.

Before running the demonstration, install all the necessary modules(RDKit, PyMol and Biopython) through creating an anaconda environment through command prompt `conda create --name myenv`  where myenv is to be replaced by own environment name. Use ` pip install chembl_webresource_client` to install ChEMBL module.

Now enter `cd -folder-name-\molecular_visualiser` to navigate to the folder of the directory where the project has been placed in.  

Then type `activate env` to open the conda environment.

Now enter `python` to open the Python shell. 

Now type the following code:

```python
>>> from protein import Protein
>>> from chem_mol import ChemicalMolecule
```

### Demonstration Workflow:
-----------

