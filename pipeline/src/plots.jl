# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Plots
using StatsPlots

include("complexity.jl")

groups = ["A", "B", "C", "D", "E", "F", "G", "H"]

grouped_datasets = Dict()

grouped_datasets["A"] = ["MR-0311", "MR-0309", "MR-0306", "MR-0303", "MR-0300-t1", "MR-0299-t2", "MR-0298-t1", "MR-0296-t2"]
grouped_datasets["B"] = ["MR-0294", "MR-0289", "MR-0288-t1", "MR-0284", "MR-0282", "MR-0276", "MR-0273", "MR-0270"]
grouped_datasets["C"] = ["MR-0313", "MR-0312", "MR-0310", "MR-0307-t2", "MR-0305", "MR-0304-t2", "MR-0302-t1", "MR-0301-t2"]
grouped_datasets["D"] = ["MR-0293", "MR-0292-t1", "MR-0291", "MR-0290", "MR-0287", "MR-0285-t1", "MR-0280-t1", "MR-0278", "MR-0275", "MR-0274"]
grouped_datasets["E"] = ["MR-0593_nd4", "MR-0592_nd4", "MR-0591_nd4", "MR-0573_nd4", "MR-0483", "MR-0460", "MR-0456"]
grouped_datasets["F"] = ["MR-0599_nd4", "MR-0597_nd4", "MR-0596_nd4", "MR-0569_nd4", "MR-0554", "MR-0465-t2"]
grouped_datasets["G"] = ["MR-0586_nd4", "MR-0585_nd4", "MR-0583_nd4", "MR-0579_nd4", "MR-0577_nd4", "MR-0575_nd4"]
grouped_datasets["H"] = ["MR-0588_nd4", "MR-0587_nd4", "MR-0582_nd4", "MR-0568_nd4", "MR-0552", "MR-0548-t1", "MR-0530"]

mse_values = Dict()
rcmse_values = Dict()

for group in groups
	mse_values[group] = []
	rcmse_values[group] = []
end

#complexity ranges
r_list = [0.2]

for group in groups
	for dataset in grouped_datasets[group]
		for r in r_list
			# open processed file
			h5open("./processed_data/"*dataset*"/"*dataset*"_processed.h5", "r") do processed_file

				plot(size=(1024, 768))
				signal = read(processed_file, "resampled_signal/data")
				plot!((0:length(signal)-1)./250, signal, label="Pipelined signal")
				plot!(xlabel="Time", ylabel="Amplitude", title=dataset)
				savefig("./processed_data/"*dataset*"/"*dataset*"_signal.png")

				plot(size=(800, 600))
				# read filtered_signal
				mse = read(processed_file, "complexity_curves/MSE_2_"*string(r))
				rcmse = read(processed_file, "complexity_curves/RCMSE_2_"*string(r))

				#compute complexity value
				mse_values[group] = [mse_values[group]; compute_complexity_value(mse)]
				rcmse_values[group] = [rcmse_values[group]; compute_complexity_value(rcmse)]

				plot!(1:80, mse, label="MSE")
				plot!(1:80, rcmse, label="RCMSE")
				plot!(xlabel="Scale", ylabel="Complexity", title=dataset)
				#y label from 0 to 1
				plot!(ylims=(0, 1))
				savefig("./processed_data/"*dataset*"/"*dataset*"_complexity_"*string(r)*".png")
			end
		end
	end
end

#violin plot
plot(size=(800, 600), legend=false)
for group in groups
	violin!([group for i in 1:length(mse_values[group])], mse_values[group], label=group, line = 0, fill = (0.2, :blue))
	boxplot!([group for i in 1:length(mse_values[group])], mse_values[group], label=group, line = (2, :black), fill = (0.3, :orange))
end
plot!(xlabel="Group", ylabel="Complexity value (80 scales)", title="Complexity distribution (MSE 2 0.2)")
savefig("./processed_data/complexity_distribution_MSE.png")

plot(size=(800, 600), legend=false)
for group in groups
	violin!([group for i in 1:length(rcmse_values[group])], rcmse_values[group], label=group, line = 0, fill = (0.2, :blue))
	boxplot!([group for i in 1:length(rcmse_values[group])], rcmse_values[group], label=group, line = (2, :black), fill = (0.3, :orange))
end
plot!(xlabel="Group", ylabel="Complexity value (80 scales)", title="Complexity distribution (RCMSE 2 0.2)")
savefig("./processed_data/complexity_distribution_RCMSE.png")
