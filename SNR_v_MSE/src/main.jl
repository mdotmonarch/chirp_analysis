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

	#=
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
	=#

	scales_list_ascendant = [
		[i for i in 1:5],
		[i for i in 1:10],
		[i for i in 1:15],
		[i for i in 1:20],
		[i for i in 1:25],
		[i for i in 1:30],
		[i for i in 1:35],
		[i for i in 1:40]
	]

	scales_list_descendant = [
		[i for i in 35:40],
		[i for i in 30:40],
		[i for i in 25:40],
		[i for i in 20:40],
		[i for i in 15:40],
		[i for i in 10:40],
		[i for i in 5:40],
		[i for i in 1:40]
	]

	scales_list_middle = [
		[i for i in 18:22],
		[i for i in 16:24],
		[i for i in 14:26],
		[i for i in 12:28],
		[i for i in 9:31],
		[i for i in 6:34],
		[i for i in 3:37],
		[i for i in 1:40]
	]

	ascendant_plots = []
	descendant_plots = []
	middle_plots = []

	h5open("./SNR_v_MSE/processed_data/"*dataset*"/"*dataset*"_processed.h5", "r") do processed_file

		f = read(processed_file)

		for scales in scales_list_ascendant
			snr_list = Float64[]
			mse_list = Float64[]
			for n in 0:251
				electrode = "electrode_"*string(n)
				snr_list = [snr_list; f[electrode]["snr"]]
				mse_list = [mse_list; compute_complexity(f[electrode]["mse"], scales)]
			end
			ascendant_plots = [ascendant_plots; scatter(snr_list, mse_list, title = dataset*": "*string(scales[1])*" to "*string(scales[end]), xlabel = "SNR", ylabel = "Complexity", size = (800, 800))]
		end

		for scales in scales_list_descendant
			snr_list = Float64[]
			mse_list = Float64[]
			for n in 0:251
				electrode = "electrode_"*string(n)
				snr_list = [snr_list; f[electrode]["snr"]]
				mse_list = [mse_list; compute_complexity(f[electrode]["mse"], scales)]
			end
			descendant_plots = [descendant_plots; scatter(snr_list, mse_list, title = dataset*": "*string(scales[1])*" to "*string(scales[end]), xlabel = "SNR", ylabel = "Complexity", size = (800, 800))]
		end

		for scales in scales_list_middle
			snr_list = Float64[]
			mse_list = Float64[]
			for n in 0:251
				electrode = "electrode_"*string(n)
				snr_list = [snr_list; f[electrode]["snr"]]
				mse_list = [mse_list; compute_complexity(f[electrode]["mse"], scales)]
			end
			middle_plots = [middle_plots; scatter(snr_list, mse_list, title = dataset*": "*string(scales[1])*" to "*string(scales[end]), xlabel = "SNR", ylabel = "Complexity", size = (800, 800))]
		end

		plot([ascendant_plots; descendant_plots; middle_plots]..., layout = (3, 8), size = (4000, 1600))
		savefig("./SNR_v_MSE/plots/"*dataset*"_snr_v_mse.png")
	end
end