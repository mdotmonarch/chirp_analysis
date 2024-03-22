# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Statistics
using Plots

include("../../lib/complexity.jl")

datasets = ["MR-0311", "MR-0294", "MR-0313", "MR-0293", "MR-0593_nd4", "MR-0599_nd4", "MR-0586_nd4", "MR-0588_nd4"]

sampling_rate = 20000
resampling_rate = 250

for dataset in datasets
	println("Processing dataset: ", dataset)

	if !isdir("./SNR_v_MSE/processed_data/"*dataset*"/")
		mkpath("./SNR_v_MSE/processed_data/"*dataset*"/")
	end

	h5open("./data/"*dataset*"_electrodes_raw.hdf5", "r") do file
		stream = read(file)
		for n in 0:251
			electrode = "electrode_"*string(n)
			signal = stream[electrode]["datos"]

			# resample signal
			step = trunc(Int, sampling_rate / resampling_rate)
			signal = signal[1:step:end]

			# calculate SNR
			noise_c = mean(signal[resampling_rate*36:resampling_rate*38].^2)
			signal_c = mean(signal.^2)
			snr = (signal_c/noise_c)

			# calculate MSE
			mse = multiscale_entropy(signal, 2, 0.2*std(signal), "sample", [i for i in 1:100])

			# open processed file and save data
			h5open("./SNR_v_MSE/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
				g = create_group(processed_file, electrode)
				g["snr"] = snr
				g["mse"] = mse
			end
		end
	end
end