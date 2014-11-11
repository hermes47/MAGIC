import os
''' characters which need to be parsed nicely '''
STRING_BREAK_CHARACTERS = "{(<"
''' Possible characters for random name generation '''
RANDOM_NAME_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
''' Length of names for MTB/PDB files '''
NAME_LENGTH = 4
''' Number of Tiers used for joining.
Caution is advised when altering this constant as OVL files will need to be updated and, 
more importantly, the use of the default 5 is fairly hard coded into MAGIC. '''
NUM_TIERS = 5
''' Root directory '''
ROOT_DIRECTORY = os.getcwd()
''' Pickled files '''
SAVED_PICKLED_FILES = os.listdir(path=ROOT_DIRECTORY+'/PickledFiles/')
''' Atomic Radii of the elements. measured in angstroms. Sourced from Wikipedia'''
ATOMIC_RADII = {"H":0.53,"C":0.67,"N":0.56,"O":0.48,"P":0.98,"F":0.42,"S":0.88,"Cl":0.79,"Cu":1.45,
                "Na":1.90, "Mg":1.45, "Si":1.11, "Ar":0.71, "Ca":1.94, "Fe":1.56, "Zn":1.42, "Br":0.94}
''' Charge redistribution rules on/off switching'''
CHARGE_RULES = {"rule1":False,"rule2":True,"rule3":False,"rule4":True,"rule5":True,"rule6":True,"rule7":True}
''' Prealign molecules'''
MOLECULE_PREALIGNMENT = True
''' Further alignment style. Choice are: dihedral, sweetspot'''
FURTHER_ALIGNMENT_STYLE = "sweetspot"

# the list of checked atoms. used for graphSearch
CHECKED_ATOMS = []

# the molecules that have been loaded into memory
LOADED_MOLECULES = {}

# Root names of all the MTB/PDB files 
MTB_PDB_ROOT_NAMES = []

