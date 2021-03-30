# coding: utf-8
import pickle
import glob

pickles = glob.glob('*.pkl')

data = {}

for file in pickles:
    with open(file, 'rb') as f:
        data[file[-10:-4]] = pickle.load(f)

with open("hedgehog.pkl", "wb") as f:
    pickle.dump(data, f)
