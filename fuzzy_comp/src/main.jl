# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

include("../../lib/complexity.jl")

using HDF5
using Plots
using Statistics

# open processed file
f = h5open("./pipeline/processed_data/MR-0530/MR-0530_processed.h5", "r")

# read data from h5 file
data = read(f, "resampled_signal/data")

#calculate mse 
mse = multiscale_entropy(data, 2, 0.2*std(data), "sample", [i for i in 1:45])

#calculate fmse
fmse = multiscale_entropy(data, 2, 0.2*std(data), "fuzzy", [i for i in 1:45])

# write to txt
open("./fuzzy_comp/processed_data/fmse_noise_julia.txt", "w") do f
	for i in axes(fmse, 1)
		write(f, string(fmse[i], "\n"))
	end
end

open("./fuzzy_comp/processed_data/mse_noise_julia.txt", "w") do f
	for i in axes(mse, 1)
		write(f, string(mse[i], "\n"))
	end
end