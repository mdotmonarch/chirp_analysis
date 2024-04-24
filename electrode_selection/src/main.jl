# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Statistics
using Plots


include("../../lib/complexity.jl")

datasets = ["MR-0596_nd4"]
sampling_rate = 20000
resampling_rate = 250

#=
for dataset in datasets
	if isdir("./electrode_selection/processed_data/"*dataset*"/")
		continue
	end
	mkpath("./electrode_selection/processed_data/"*dataset*"/")
	local f = h5open("./data/"*dataset*"_electrodes_raw.hdf5", "r")
	stream = read(f)
	close(f)

	for n in 0:251
		println("Processing electrode: ", n)
		electrode = "electrode_"*string(n)
		signal = stream[electrode]["datos"]

		#resample signal
		step = trunc(Int, sampling_rate / resampling_rate)
		signal = signal[1:step:end]

		# calculate SNR
		noise_c = mean(signal[resampling_rate*36:resampling_rate*38].^2)
		signal_c = mean(signal.^2)
		snr = (signal_c/noise_c)

		# calculate fRCMSE
		fRCMSE = refined_composite_multiscale_entropy(signal, 2, 0.2*std(signal), "fuzzy", [i for i in 20:45])
		cmp = compute_complexity(fRCMSE)

		# open processed file and save data
		p_f = h5open("./electrode_selection/processed_data/"*dataset*"/"*dataset*"_processed.h5", "cw")
		g = create_group(p_f, electrode)
		g["snr"] = snr
		g["cmp"] = cmp
		close(p_f)
	end
end
=#
for dataset in datasets
	p_f = h5open("./electrode_selection/processed_data/"*dataset*"/"*dataset*"_processed.h5", "r")
	data = read(p_f)
	close(p_f)

	snr = Float64[]
	cmp = Float64[]

	for n in 0:251
		electrode = "electrode_"*string(n)
		snr = [snr; data[electrode]["snr"]]
		cmp = [cmp; data[electrode]["cmp"]]
	end

	plot(size=(800, 600), xlims=(0, maximum(snr)+5), ylims=(0, maximum(cmp)+0.1))
	scatter!(snr, cmp, label="Electrodes", xlabel="SNR", ylabel="CMP", title=dataset*": Electrode selection", legend=:topright, legendfontsize=8, legendtitlefontsize=8, grid=true, framestyle = :box, xticks=0:5:maximum(snr), yticks=0:0.1:maximum(cmp))
	savefig("./electrode_selection/plots/"*dataset*"_electrode_distribution.png")

	reg_data = [snr cmp]
	reg_data = reg_data[sortperm(reg_data[:,1]), :] # sorts data by SNR

	s_mem = Float64[]
	for i in axes(reg_data, 1)
		if i == 1
			s_mem = [s_mem; 0]
		else
			s_mem = [s_mem; s_mem[i-1] + (0.5)*(reg_data[i, 2]+reg_data[i-1, 2])*(reg_data[i, 1]-reg_data[i-1, 1])]
		end
	end

	# regression w/ function a + b*exp(c*x)

	sum_x_2 = sum((reg_data[:, 1].-reg_data[1, 1]).^2)
	sum_x_s = sum((reg_data[:, 1].-reg_data[1, 1]).*s_mem[1:end])
	sum_s_2 = sum(s_mem[1:end].^2)
	sum_y_x = sum((reg_data[:, 2].-reg_data[1, 2]).*(reg_data[:, 1].-reg_data[1, 1]))
	sum_y_s = sum((reg_data[:, 2].-reg_data[1, 2]).*s_mem[1:end])

	M_1 = [sum_x_2 sum_x_s; sum_x_s sum_s_2]
	Y_1 = [sum_y_x; sum_y_s]

	AB = inv(M_1)*Y_1

	A = AB[1]
	B = AB[2]

	sum_theta = sum(exp.(B*(reg_data[:, 1])))
	sum_theta_2 = sum(exp.(2*B*(reg_data[:, 1])))
	sum_y = sum(reg_data[:, 2])
	sum_y_theta = sum(reg_data[:, 2].*exp.(B*reg_data[:, 1]))

	M_2 = [length(snr) sum_theta; sum_theta sum_theta_2]
	Y_2 = [sum_y; sum_y_theta]

	ab = inv(M_2)*Y_2

	a = ab[1]
	b = ab[2]
	c = B

	println("Regression approximation: ")
	println("a: ", a)
	println("b: ", b)
	println("c: ", c)
	println("")

	# find elbow point
	
	x = [i for i in range(0, stop=maximum(snr), length=1000)]
	y = (b*exp.(c*(x))).+a

	x_n = (x[:].-minimum(x))./(maximum(x)-minimum(x))
	y_n = (y[:].-minimum(y))./(maximum(y)-minimum(y))

	Dd_p = [x_n y_n.+x_n]
	Dd_m = [x_n y_n.-x_n]

	plot(size=(800, 600))
	plot!(Dd_p[:, 1], Dd_p[:, 2])
	savefig("./electrode_selection/plots/"*dataset*"_data_elbow_p.png")

	plot(size=(800, 600))
	plot!(Dd_m[:, 1], Dd_m[:, 2])
	savefig("./electrode_selection/plots/"*dataset*"_data_elbow_m.png")

	elbow_index = Int64

	for i in axes(Dd_p, 1)
		if i == 1 || i == size(Dd_p, 1)
			continue
		end
		if b > 0 && Dd_p[i, 2] < Dd_p[i-1, 2] && Dd_p[i, 2] < Dd_p[i+1, 2]
			println("SNR threshold: ", x[i])
			elbow_index = i
		elseif b < 0 && Dd_p[i, 2] > Dd_p[i-1, 2] && Dd_p[i, 2] > Dd_p[i+1, 2]
			println("SNR threshold: ", x[i])
			elbow_index = i
		end
	end

	for i in axes(Dd_m, 1)
		if i == 1 || i == size(Dd_m, 1)
			continue
		end
		if b > 0 && Dd_m[i, 2] < Dd_m[i-1, 2] && Dd_m[i, 2] < Dd_m[i+1, 2]
			println("SNR threshold: ", x[i])
			elbow_index = i
		elseif b < 0 && Dd_m[i, 2] > Dd_m[i-1, 2] && Dd_m[i, 2] > Dd_m[i+1, 2]
			println("SNR threshold: ", x[i])
			elbow_index = i
		end
	end

	# exponential regression
	plot(size=(800, 600), xlims=(0, maximum(snr)+5), ylims=(0, maximum(cmp)+0.1))
	scatter!(snr, cmp, label="Electrodes", xlabel="SNR", ylabel="CMP", title=dataset*": Electrode selection", legend=:topright, legendfontsize=8, legendtitlefontsize=8, grid=true, framestyle = :box, xticks=0:5:maximum(snr), yticks=0:0.1:maximum(cmp))
	plot!(x, y, label="Exponential regression curve")
	vline!([x[elbow_index]], label="SNR threshold: "*string(round(x[elbow_index], digits=3)))
	savefig("./electrode_selection/plots/"*dataset*"_electrode_selection.png")
end