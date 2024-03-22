# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Statistics
using Plots

include("../../lib/complexity.jl")

datasets = ["MR-0311", "MR-0294", "MR-0313", "MR-0293"]

for dataset in datasets
	println("Processing dataset: ", dataset)

	h5open("./SNR_v_MSE/processed_data/"*dataset*"/"*dataset*"_processed.h5", "r") do file
		stream = read(file)
		println(length(keys(stream)))
	end
end