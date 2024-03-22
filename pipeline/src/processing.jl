# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

include("../../lib/complexity.jl")

using HDF5
using DSP
using JSON
using Statistics
using Plots

function high_pass_filter(signal, nco, order)
	filter = digitalfilter(Highpass(nco), Butterworth(order))
	return filtfilt(filter, signal)
end

function band_pass_filter(signal, nco_low, nco_high, order)
	filter = digitalfilter(Bandpass(nco_low, nco_high), Butterworth(order))
	return filtfilt(filter, signal)
end

# parameters
sampling_rate = 20000
nyquist_frequency = 0.5 * sampling_rate

function select_with_signal_to_noise_ratio(dataset, snr_threshold = 10)
	# open processed file
	h5open("./analysis/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if clean_electrodes group exists
		if "clean_electrodes" in keys(read(processed_file))
			println("Skipping signal to noise ratio selection.")
			return
		end

		println("Selecting electrodes with signal to noise ratio > "*string(snr_threshold)*"... ")
		clean_electrodes = String[]
		h5open("./data/"*dataset*"_electrodes_raw.hdf5", "r") do file
			stream = read(file)
			for n in 0:251
				electrode = "electrode_"*string(n)
				signal = stream[electrode]["datos"]
				noise_c = mean(signal[sampling_rate*36:sampling_rate*38].^2)
				signal_c = mean(signal.^2)
		
				snr = (signal_c/noise_c)
				if snr > snr_threshold
					println("SNR of "*electrode*": "*string(snr))
					clean_electrodes = [clean_electrodes; electrode]
				end
			end
		end

		if length(clean_electrodes) == 0
			println("No electrodes with signal to noise ratio > "*string(snr_threshold)*".")
			return
		end

		g = create_group(processed_file, "clean_electrodes")
		g["data"] = clean_electrodes

		sg = create_group(g, "meta")
		sg["snr_threshold"] = snr_threshold
		println("Done.")
	end
end

function normalize_signals_and_average(dataset)
	# open processed file
	h5open("./analysis/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if average_signal group exists
		if "average_signal" in keys(read(processed_file))
			println("Skipping signal averaging.")
			return
		end
		# check if it has clean electrodes
		if !("clean_electrodes" in keys(read(processed_file)))
			return
		end


		println("Normalizing selected signals and averaging... ")
		# read clean electrodes data
		clean_electrodes = read(processed_file, "clean_electrodes/data")
		# read raw data
		h5open("./data/"*dataset*"_electrodes_raw.hdf5", "r") do file
			stream = read(file)
			# create normalized signals group
			for electrode in clean_electrodes
				# normalize signal and store it in the processed file
				signal = stream[electrode]["datos"]
				max_signal = maximum(signal)
				min_signal = minimum(signal)
				normalized_signal = ((signal .- min_signal) ./ (max_signal - min_signal))
				normalized_signal = (2 .* normalized_signal) .- 1

				if electrode == clean_electrodes[1]
					average_signal = normalized_signal
				else
					average_signal = average_signal + normalized_signal
				end
			end
			average_signal = average_signal ./ length(clean_electrodes)
			g = create_group(processed_file, "average_signal")
			g["data"] = average_signal
			println("Done.")
		end
	end
end

function filter_signal(dataset, filter)
	# open processed file
	h5open("./analysis/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if average_signal group exists
		if "filtered_signal" in keys(read(processed_file))
			println("Skipping signal filtering.")
			return
		end
		# check if it has clean electrodes
		if !("clean_electrodes" in keys(read(processed_file)))
			return
		end

		println("Filtering signal... ")
		# read average_signal
		average_signal = read(processed_file, "average_signal/data")

		# filter signal
		if filter["type"] == "HPF"
			filtered_signal = high_pass_filter(average_signal, filter["cutoff"]/nyquist_frequency, 5)
		elseif filter["type"] == "BPF"
			filtered_signal = band_pass_filter(average_signal, filter["cutoff_low"]/nyquist_frequency, filter["cutoff_high"]/nyquist_frequency, 5)
		end

		g = create_group(processed_file, "filtered_signal")
		g["data"] = filtered_signal
		
		sg = create_group(g, "meta")
		sg["filter_type"] = filter["type"]
		if filter["type"] == "HPF"
			sg["cutoff"] = filter["cutoff"]
		elseif filter["type"] == "BPF"
			sg["cutoff_low"] = filter["cutoff_low"]
			sg["cutoff_high"] = filter["cutoff_high"]
		end
		println("Done.")
	end
end

function resample_signal(dataset, resampling_rate)
	# open processed file
	h5open("./analysis/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if average_signal group exists
		if "resampled_signal" in keys(read(processed_file))
			println("Skipping signal resampling.")
			return
		end
		# check if it has clean electrodes
		if !("clean_electrodes" in keys(read(processed_file)))
			return
		end

		println("Resampling signal... ")
		# read filtered_signal
		filtered_signal = read(processed_file, "filtered_signal/data")

		# resample signal
		step = trunc(Int, sampling_rate / resampling_rate)
		resampled_signal = filtered_signal[1:step:end]

		g = create_group(processed_file, "resampled_signal")
		g["data"] = resampled_signal

		sg = create_group(g, "meta")
		sg["resampling_rate"] = resampling_rate
		println("Done.")
	end
end

function compute_complexity_curve(dataset, type, m, r, scales)
	# open processed file
	h5open("./analysis/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if it has clean electrodes
		if !("clean_electrodes" in keys(read(processed_file)))
			return
		end

		println("Computing "*type*" complexity curve... ")
		# read resampled_signal
		resampled_signal = read(processed_file, "resampled_signal/data")

		# compute complexity curve
		if type == "MSE"
			complexity_curve = multiscale_entropy(resampled_signal, m, r, "sample", scales)
		elseif type == "RCMSE"
			complexity_curve = refined_composite_multiscale_entropy(resampled_signal, m, r, "sample", scales)
		elseif type == "FMSE"
			complexity_curve = multiscale_entropy(resampled_signal, m, r, "fuzzy", scales)
		elseif type == "FRCMSE"
			complexity_curve = refined_composite_multiscale_entropy(resampled_signal, m, r, "fuzzy", scales)
		end

		g = create_group(processed_file, type)
		g["data"] = complexity_curve

		sg = create_group(g, "meta")
		sg["m"] = m
		sg["r"] = r
		sg["scales"] = scales
		println("Done.")
	end
end