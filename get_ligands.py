import os
import csv
import pandas as pd
from shutil import copyfile


with open('info.csv') as csv_file:
    pdb_info = csv.reader(csv_file, delimiter=',')
    all_li = []
    for line in pdb_info:
        pdb = line[0].lower()
        lig_id = line[3].split('_')[-1].lower()
        if line[2] == 't':
            src = "/home/cradoux/Research/Hotspots/validation2/validation/ligands/{}.mol2".format(pdb)
            dst = "pdb_files/{}/{}/ligand.mol2".format(pdb[1:3], pdb)
            try:
                copyfile(src=src,dst=dst)
            except IOError:
                continue