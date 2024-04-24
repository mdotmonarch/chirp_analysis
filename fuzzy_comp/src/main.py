# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

import h5py
import neurokit2 as nk

#open processed data
f = h5py.File("./pipeline/processed_data/MR-0530/MR-0530_processed.h5", "r")

#read data from h5 file
data = f['resampled_signal/data']

#calculate mse 
mse = nk.complexity_mse(data, scale=[i for i in range(1, 46)], dimension=2)
mse = mse[1]["Value"]
print()

#calculate fmse
fmse = nk.complexity_mse(data, scale=[i for i in range(1, 46)], dimension=2, fuzzy=True)
fmse = fmse[1]["Value"]

#write data to text
with open("./fuzzy_comp/processed_data/mse_noise_python.txt", "w") as file:
	for i in range(len(mse)):
		file.write(str(mse[i]) + "\n")

with open("./fuzzy_comp/processed_data/fmse_noise_python.txt", "w") as file:
	for i in range(len(fmse)):
		file.write(str(fmse[i]) + "\n")