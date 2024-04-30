# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

include("../../lib/complexity.jl")

using HDF5
using DSP
using JSON
using Statistics
using Plots
using CurveFit

function group_check(file, list, i=0)
	if i == length(list)
		return true
	end
	path = "/"*join(list[1:i], "/")
	println("Checking path: ", path)
	if list[i+1] in keys(file[path])
		return group_check(file, list, i+1)
	end
	return false
end

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
resampling_rate = 250
nyquist_frequency = 0.5 * sampling_rate

function use_all_electrodes(dataset)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if clean_electrodes group exists
		if group_check(processed_file, ["clean_electrodes"])
			println("Skipping clean electrodes assignment...")
			return
		end

		println("Using all electrodes...")
		clean_electrodes = String[]

		for n in 0:251
			electrode = "electrode_"*string(n)
			clean_electrodes = [clean_electrodes; electrode]
		end

		processed_file["clean_electrodes/data"] = clean_electrodes
		println("Done.")
	end
end

function normalize_signals_and_average(dataset)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if average_signal group exists
		if group_check(processed_file, ["average_signal"])
			println("Skipping signal normalization and averaging.")
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
				# normalize signal using z score
				signal = stream[electrode]["datos"]
				
				u = mean(signal)
				s = std(signal)
				normalized_signal = (signal .- u) ./ s

				if electrode == clean_electrodes[1]
					average_signal = normalized_signal
				else
					average_signal = average_signal + normalized_signal
				end
			end
			average_signal = average_signal ./ length(clean_electrodes)

			processed_file["average_signal/data"] = average_signal
			println("Done.")
		end
	end
end

function filter_signal(dataset, filter)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if average_signal group exists
		if group_check(processed_file, ["filtered_signal"])
			println("Skipping signal filtering.")
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

		processed_file["filtered_signal/data"] = filtered_signal
		processed_file["filtered_signal/meta/filter_type"] = filter["type"]
		if filter["type"] == "HPF"
			processed_file["filtered_signal/meta/cutoff"] = filter["cutoff"]
		elseif filter["type"] == "BPF"
			processed_file["filtered_signal/meta/cutoff_low"] = filter["cutoff_low"]
			processed_file["filtered_signal/meta/cutoff_high"] = filter["cutoff_high"]
		end
		println("Done.")
	end
end

function resample_signal(dataset, resampling_rate)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if average_signal group exists
		if group_check(processed_file, ["signals", "complete"])
			println("Skipping signal to noise ratio selection.")
			return
		end

		println("Resampling signal... ")
		# read filtered_signal
		filtered_signal = read(processed_file, "filtered_signal/data")

		# resample signal
		step = trunc(Int, sampling_rate / resampling_rate)
		resampled_signal = filtered_signal[1:step:end]

		processed_file["signals/complete/data"] = resampled_signal
		processed_file["signals/meta/resampling_rate"] = resampling_rate
		println("Done.")
	end
end

function get_signal_segments(dataset)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if segments exist
		if group_check(processed_file, ["signals", "flash"])
			println("Skipping flash segment extraction.")
			return
		else
			signal = read(processed_file, "signals/complete/data")
			processed_file["signals/flash/data"] = signal[1:6*250]
		end

		if group_check(processed_file, ["signals", "frequency"])
			println("Skipping increasing frequency segment extraction.")
			return
		else
			signal = read(processed_file, "signals/complete/data")
			processed_file["signals/frequency/data"] = signal[8*250:23*250]
		end

		if group_check(processed_file, ["signals", "amplitude"])
			println("Skipping increasing amplitude segment extraction.")
			return
		else
			signal = read(processed_file, "signals/complete/data")
			processed_file["signals/amplitude/data"] = signal[25*250:33*250]
		end
	end
end

function compute_complexity_curve(dataset, segment, type, m, r, scales)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if complexity group exists
		if group_check(processed_file, ["complexity", segment, type, string(r)])
			println("Skipping signal to noise ratio selection.")
			return
		end

		println("Computing complexity curve... ")

		# read signal
		signal = read(processed_file, "signals/"*segment*"/data")

		# compute complexity curve
		if type == "MSE"
			complexity_curve = multiscale_entropy(signal, m, r*std(signal), "sample", scales)
		elseif type == "RCMSE"
			complexity_curve = refined_composite_multiscale_entropy(signal, m, r*std(signal), "sample", scales)
		elseif type == "FMSE"
			complexity_curve = multiscale_entropy(signal, m, r*std(signal), "fuzzy", scales)
		elseif type == "FRCMSE"
			complexity_curve = refined_composite_multiscale_entropy(signal, m, r*std(signal), "fuzzy", scales)
		end

		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/data"] = complexity_curve
		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/m"] = m
		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/r"] = r
		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/scales"] = scales
	end
end

function compute_linear_regression(dataset, segment, type, r)
	# open processed file
	h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw") do processed_file
		# check if complexity group exists
		if group_check(processed_file, ["complexity", segment, type, string(r), "meta", "linear_fit"])
			println("Skipping linear regression.")
			return
		end

		println("Computing linear regression... ")

		# read complexity curve
		complexity_curve = read(processed_file, "complexity/"*segment*"/"*type*"/"*string(r)*"/data")
		scales = read(processed_file, "complexity/"*segment*"/"*type*"/"*string(r)*"/meta/scales")

		# compute linear regression
		a, b = linear_fit(scales, complexity_curve)

		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/linear_fit/a"] = a
		processed_file["complexity/"*segment*"/"*type*"/"*string(r)*"/meta/linear_fit/b"] = b
	end
end