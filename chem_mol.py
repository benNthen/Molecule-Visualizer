from __future__ import print_function

import pandas as pd
import os
import rdkit
import sys
import pymol

from os import path
from rdkit import Chem
from pymol import cmd
from pathlib import Path
from rdkit.Chem import AllChem
from chembl_webresource_client.new_client import new_client

substructure = new_client.substructure


# A ChemicalMolecule class containing all the necessary small chemical molecule information for the visualiser
class ChemicalMolecule:

    def __init__(self, chemblid):
        self.id = chemblid

        # Search ChEMBLdb for molecule matching the chem id provided
        # then retrieves its molecular structure
        res = substructure.filter(chembl_id=chemblid).only(['molecule_structures'])

        # Store the molecular structure into a dataframe
        # then extract only the canonical smiles attribute
        mol_df = pd.DataFrame(res).values[0][0]
        smiles = list(mol_df.values())[0]
        mol = Chem.MolFromSmiles(smiles)

        # Generate 3D co-ordinates based on canonical smiles
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)

        # save the 3D co-ordinates into a .mol file
        filename = chemblid + '.mol'

        # set up destination of created mol file to be inside the data folder of the project's directory
        current_path = os.getcwd() + '/' + filename
        folder_path = os.getcwd() + r'/data/' + filename

        # Before creating the mol file, the data folder is checked if a similar file already exists
        # If not, then the newly created file is immediately moved to the data folder
        file_exist = path.isfile(folder_path)

        if not file_exist:
            rdkit.Chem.rdmolfiles.MolToMolFile(mol, filename)

            os.rename(current_path, folder_path)

        # Displays a message if data retrieval and mol file creation was successful
        self.my_file = Path(os.getcwd() + r'/data/' + filename)

        if self.my_file.is_file():
            print('A ChEMBL matching ' + chemblid + ' was found and saved in ' + str(self.my_file))

        else:
            print('No ChEMBL matching' + chemblid + ' was found. Please try again.')

    # a show instance method which will display the molecule on request in the running demonstration via PyMol
    def show(self):

        # Ensures that show instance will only work if the model exists in the directory
        # otherwise display a message prompt error
        if self.my_file.is_file():

            print('Processing molecule model for ' + self.id + ' to PyMol')

            # Launches PyMol via Python command script
            _stdouterr = sys.stdout, sys.stderr
            pymol.finish_launching(['pymol', '-q'])
            sys.stdout, sys.stderr = _stdouterr

            # Load the mol file into PyMol to view its 3D model
            cmd.get_view(0)
            cmd.load(self.my_file)

        else:
            print('Error. No model was loaded to PyMol.')

        return self.id
