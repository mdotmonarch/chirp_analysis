# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Plots
using StatsPlots
using FFTW

include("../../lib/complexity.jl")

function signal_reconstruct(signal, info)
	return [(s-info[:ADZero])*(info[:ConversionFactor]*(10.0^info[:Exponent])) for s in signal]
end

#=

file = h5open("noise_analysis/data/chirp.h5", "r")
stream = read(file, "Data/Recording_0/AnalogStream/Stream_0")
close(file)

h5open("./noise_analysis/processed_data/noise_chirp.h5", "cw") do processed_file
	signal_length = trunc(Int, length(stream["ChannelData"])/252)	# 252 channels
	for i in 1:252
		electrode_label = "electrode_"*stream["InfoChannel"][i][:Label]
		println("Processing electrode #"*string(i-1)*": "*electrode_label)

		info = stream["InfoChannel"][i]
		signal_raw = stream["ChannelData"][signal_length*(i-1)+1:signal_length*i]

		signal = signal_reconstruct(signal_raw, info)
		signal = [s*1000 for s in signal] # Convert from V to mV

		processed_file["electrode_"*string(i-1)*"/datos"] = signal
	end
end
close(file)
=#

h5open("./noise_analysis/processed_data/noise_chirp.h5", "r") do processed_file
	for i in 0:251
		electrode_label = "electrode_"*string(i)*"/datos"
		signal = read(processed_file, electrode_label)

		#resample signal
		step = trunc(Int, 20000 / 250)
		resampled_signal = signal[1:step:end]
		time_plot = plot(resampled_signal, title = "Signal", label='f', legend=:top)

		complexity = 

		#compute fft
		F = fftshift(fft(signal))
		freqs = fftshift(fftfreq(length(signal), 20000))
		freq_plot = plot(freqs, abs.(F), title = "Spectrum", xlim=(-10000, +10000), xticks=-10000:2000:10000, label="abs.(F)", legend=:top)

		plot(time_plot, freq_plot, layout = (2,1), size=(4000, 2000))
		savefig("./noise_analysis/plots/electrode_"*string(i)*".png")
	end
end