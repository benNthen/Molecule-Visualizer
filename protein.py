import os
import sys
import pymol

from pymol import cmd
from pathlib import Path
from Bio.PDB import PDBList


# A Protein class containing all the necessary protein information for the visualiser
class Protein:

    def __init__(self, pdbid):

        self.id = pdbid

        pdbl = PDBList()

        # Extract project folder's current directory address
        dir = os.getcwd()
        folder_path = dir + r'\data'

        # Search the Protein Data Bank for matching ID, if model exists then
        # download and save on the current directory as mmCif file
        protein = pdbl.retrieve_pdb_file(pdbid, file_format='mmCif', pdir=folder_path)

        filename = pdbid + '.cif'

        self.my_file = Path(folder_path + '/' + filename)

        # Displays a message if data retrieval and mol file creation was successful
        if self.my_file.is_file():
            print('A protein with an ID matching ' + pdbid + ' was found and saved in ' + str(self.my_file))

        else:
            print('Protein ID entered is invalid. Please try again.')

    # A show instance method which will display the molecule on request in the running demonstration
    def show(self):

        # Ensures that show instance will only work if the model exists in the directory
        # otherwise display a message prompt error
        if self.my_file.is_file():

            print('Processing protein model for ' + self.id + ' to PyMol')

            # Launches PyMol via Python command script
            _stdouterr = sys.stdout, sys.stderr
            pymol.finish_launching(['pymol', '-q'])
            sys.stdout, sys.stderr = _stdouterr

            # Load the mmCif file on PyMol to view its 3D model
            cmd.get_view(0)
            cmd.load(self.my_file)

        else:
            print('Error. No model was loaded to PyMol.')

        return self.id
