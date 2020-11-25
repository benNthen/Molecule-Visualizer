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
The following will demonstrate how the workflow will retrieve molecular data online and convert to be visualised as a 3D object on PyMol simulatenously. 

```python
# Get the crystal structure of the main COVID-19 protease
>>> protease = Protein(pdbid='6LU7')
Downloading PDB structure '6LU7'...
A protein with an ID matching 6LU7 was found and saved in C:\folder-name\molecular_visualiser\data\6lu7.cif
>>>protease.show()
Processing protein model for 6LU7 to PyMol
 ExecutiveLoad-Detail: Detected mmCIF
'6LU7'
>>>  Detected OpenGL version 4.6. Shaders available.
 Detected GLSL version 4.60.
 License Expiry date: 01-jun-2021
```
[figure 1 image here]

As shown in Figure 1, this workflow is able to retrieve data online and visual molecular structures in 3D using PyMol. 

```python
# Get the 3D structure of lopinavir
>>> lopinavir = ChemicalMolecule(chemblid='CHEMBL729')
A ChEMBL matching CHEMBL729 was found and saved in C:\Users\bened\Desktop\molecular_visualiser\data\CHEMBL729.mol
>>>lopinavir.show()     
Processing molecule model for CHEMBL729 to PyMol
'CHEMBL729'
```
[insert figure 2 here]

Figure 2 shows that it is possible to add a second molecular object to the same PyMol window simulatenously. 
