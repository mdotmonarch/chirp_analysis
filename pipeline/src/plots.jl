# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Plots
using StatsPlots

include("../../lib/complexity.jl")

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
fmse_values = Dict()
frcmse_values = Dict()

for group in groups
	mse_values[group] = []
	rcmse_values[group] = []
	fmse_values[group] = []
	frcmse_values[group] = []
end

for group in groups
	for dataset in grouped_datasets[group]
		println("Plotting dataset: ", dataset)
		# if directory does not exist, create it
		if !isdir("./pipeline/plots/"*dataset*"/")
			mkpath("./pipeline/plots/"*dataset*"/")
		end

		# open processed file
		h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "r") do processed_file
			if !("resampled_signal" in keys(read(processed_file)))
				return
			end

			plot(size=(1024, 768))
			signal = read(processed_file, "resampled_signal/data")
			plot!((0:length(signal)-1)./250, signal, label="Pipelined signal")
			plot!(xlabel="Time", ylabel="Amplitude", title=dataset)
			plot!(ylims=(-1, 1))
			savefig("./pipeline/plots/"*dataset*"/"*dataset*"_signal.png")

			plot(size=(800, 600))
			# read entropy curves
			mse_curve = read(processed_file, "MSE/data")
			mse_scales = read(processed_file, "MSE/meta/scales")
			rcmse_curve = read(processed_file, "RCMSE/data")
			rcmse_scales = read(processed_file, "RCMSE/meta/scales")
			fmse_curve = read(processed_file, "FMSE/data")
			fmse_scales = read(processed_file, "FMSE/meta/scales")
			frcmse_curve = read(processed_file, "FRCMSE/data")
			frcmse_scales = read(processed_file, "FRCMSE/meta/scales")

			#compute complexity value
			mse_values[group] = [mse_values[group]; compute_complexity(mse_curve, mse_scales)]
			rcmse_values[group] = [rcmse_values[group]; compute_complexity(rcmse_curve, rcmse_scales)]
			fmse_values[group] = [fmse_values[group]; compute_complexity(fmse_curve, fmse_scales)]
			frcmse_values[group] = [frcmse_values[group]; compute_complexity(frcmse_curve, frcmse_scales)]

			plot!(mse_scales, mse_curve, label="MSE")
			plot!(rcmse_scales, rcmse_curve, label="RCMSE")
			plot!(fmse_scales, fmse_curve, label="FMSE")
			plot!(frcmse_scales, frcmse_curve, label="FRCMSE")
			
			plot!(xlabel="Scale", ylabel="Entropy", title=dataset)
			#y label from 0 to 1
			plot!(ylims=(0, 0.2))
			savefig("./pipeline/plots/"*dataset*"/"*dataset*"_entropy.png")
		end
	end
end

println("Plotting complexity distribution...")
println("# of group A datasets: ", length(mse_values["A"]))
println("# of group B datasets: ", length(mse_values["B"]))
println("# of group C datasets: ", length(mse_values["C"]))
println("# of group D datasets: ", length(mse_values["D"]))
println("# of group E datasets: ", length(mse_values["E"]))
println("# of group F datasets: ", length(mse_values["F"]))
println("# of group G datasets: ", length(mse_values["G"]))
println("# of group H datasets: ", length(mse_values["H"]))

#violin plot
plot(size=(800, 600), legend=false)
for group in groups
	violin!([group for i in 1:length(mse_values[group])], mse_values[group], label=group, line = 0, fill = (0.2, :blue))
	boxplot!([group for i in 1:length(mse_values[group])], mse_values[group], label=group, line = (2, :black), fill = (0.3, :orange))
end
plot!(xlabel="Group", ylabel="Complexity", title="Complexity distribution (MSE)")
savefig("./pipeline/plots/complexity_distribution_MSE.png")

plot(size=(800, 600), legend=false)
for group in groups
	violin!([group for i in 1:length(rcmse_values[group])], rcmse_values[group], label=group, line = 0, fill = (0.2, :blue))
	boxplot!([group for i in 1:length(rcmse_values[group])], rcmse_values[group], label=group, line = (2, :black), fill = (0.3, :orange))
end
plot!(xlabel="Group", ylabel="Complexity", title="Complexity distribution (RCMSE)")
savefig("./pipeline/plots/complexity_distribution_RCMSE.png")

plot(size=(800, 600), legend=false)
for group in groups
	violin!([group for i in 1:length(fmse_values[group])], fmse_values[group], label=group, line = 0, fill = (0.2, :blue))
	boxplot!([group for i in 1:length(fmse_values[group])], fmse_values[group], label=group, line = (2, :black), fill = (0.3, :orange))
end
plot!(xlabel="Group", ylabel="Complexity", title="Complexity distribution (FMSE)")
savefig("./pipeline/plots/complexity_distribution_FMSE.png")

plot(size=(800, 600), legend=false)
for group in groups
	violin!([group for i in 1:length(frcmse_values[group])], frcmse_values[group], label=group, line = 0, fill = (0.2, :blue))
	boxplot!([group for i in 1:length(frcmse_values[group])], frcmse_values[group], label=group, line = (2, :black), fill = (0.3, :orange))
end
plot!(xlabel="Group", ylabel="Complexity", title="Complexity distribution (FRCMSE)")
savefig("./pipeline/plots/complexity_distribution_FRCMSE.png")