# coding: utf-8
###################################################
# Repackage the individual pickle files produced by
# sightline_analysis(_subproc_manager).py for each
# dataset into one pickle file
###################################################

import pickle
import glob

pickles = glob.glob('*mass.pkl')

data = {}

for file in pickles:
    with open(file, 'rb') as f:
        data[file[-10:-4]] = pickle.load(f)

with open("hedgehog_mass.pkl", "wb") as f:
    pickle.dump(data, f)
