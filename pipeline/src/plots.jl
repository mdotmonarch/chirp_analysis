# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using HDF5
using Plots
using StatsPlots

include("../../lib/complexity.jl")

g_list = ["A", "B", "C", "D", "E", "F", "G", "H"]
r_list = ["0.1", "0.2", "0.3", "0.4", "0.5"]
s_list = ["complete", "flash", "frequency", "amplitude"]

group_labels = Dict()
group_labels["A"] = "WT young"
group_labels["B"] = "WT adult"
group_labels["C"] = "5xFAD young"
group_labels["D"] = "5xFAD adult"
group_labels["E"] = "XBP1s young"
group_labels["F"] = "XBP1s adult"
group_labels["G"] = "Double young"
group_labels["H"] = "Double adult"
group_labels_plot = ["WT young" "WT adult" "5xFAD young" "5xFAD adult" "XBP1s young" "XBP1s adult" "Double young" "Double adult"]

grouped_datasets = Dict()
grouped_datasets["A"] = ["MR-0311", "MR-0309", "MR-0306", "MR-0303", "MR-0300-t1", "MR-0299-t2", "MR-0298-t1", "MR-0296-t2"]
grouped_datasets["B"] = ["MR-0294", "MR-0289", "MR-0288-t1", "MR-0284", "MR-0283-t1", "MR-0282", "MR-0276", "MR-0273", "MR-0270"]
grouped_datasets["C"] = ["MR-0313", "MR-0312", "MR-0310", "MR-0307-t2", "MR-0305", "MR-0304-t2", "MR-0302-t1", "MR-0301-t2", "MR-0297"]
grouped_datasets["D"] = ["MR-0293", "MR-0292-t1", "MR-0291", "MR-0290", "MR-0287", "MR-0285-t1", "MR-0280-t1", "MR-0278", "MR-0275", "MR-0274"]
grouped_datasets["E"] = ["MR-0593_nd4", "MR-0592_nd4", "MR-0591_nd4", "MR-0573_nd4", "MR-0483", "MR-0460", "MR-0456"]
grouped_datasets["F"] = ["MR-0599_nd4", "MR-0597_nd4", "MR-0596_nd4", "MR-0569_nd4", "MR-0554", "MR-0465-t2"]
grouped_datasets["G"] = ["MR-0586_nd4", "MR-0585_nd4", "MR-0583_nd4", "MR-0579_nd4", "MR-0577_nd4", "MR-0575_nd4"]
grouped_datasets["H"] = ["MR-0588_nd4", "MR-0587_nd4", "MR-0582_nd4", "MR-0568_nd4", "MR-0552", "MR-0548-t1", "MR-0530"]

sex = Dict()
sex["M"] = [
	"MR-0311",
	"MR-0309",
	"MR-0306",
	"MR-0299-t2",
	"MR-0298-t1",
	"MR-0296-t2",
	"MR-0282",
	"MR-0276",
	"MR-0273",
	"MR-0270",
	"MR-0310",
	"MR-0305",
	"MR-0293",
	"MR-0292-t1",
	"MR-0280-t1",
	"MR-0278",
	"MR-0275",
	"MR-0592_nd4",
	"MR-0591_nd4",
	"MR-0483",
	"MR-0460",
	"MR-0456",
	"MR-0599_nd4",
	"MR-0597_nd4",
	"MR-0596_nd4",
	"MR-0569_nd4",
	"MR-0554",
	"MR-0465-t2",
	"MR-0575_nd4",
	"MR-0582_nd4",
	"MR-0552"
]
sex["F"] = [
	"MR-0303",
	"MR-0300-t1",
	"MR-0294",
	"MR-0289",
	"MR-0288-t1",
	"MR-0284",
	"MR-0283-t1",
	"MR-0313",
	"MR-0312",
	"MR-0307-t2",
	"MR-0304-t2",
	"MR-0302-t1",
	"MR-0301-t2",
	"MR-0297",
	"MR-0291",
	"MR-0290",
	"MR-0287",
	"MR-0285-t1",
	"MR-0274",
	"MR-0593_nd4",
	"MR-0573_nd4",
	"MR-0586_nd4",
	"MR-0585_nd4",
	"MR-0583_nd4",
	"MR-0579_nd4",
	"MR-0577_nd4",
	"MR-0588_nd4",
	"MR-0587_nd4",
	"MR-0568_nd4",
	"MR-0548-t1",
	"MR-0530"
]


for group in g_list
	for dataset in grouped_datasets[group]
		println("Plotting dataset: ", dataset)
		# if directory does not exist, create it
		if !isdir("./pipeline/plots/"*dataset*"/")
			mkpath("./pipeline/plots/"*dataset*"/")
		end

		# open processed file
		h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "r") do processed_file

			for segment in s_list
				plot(size=(1024, 768))
				signal = read(processed_file, "signals/"*segment*"/data")
				plot!((0:length(signal)-1)./250, signal, label="Response")
				plot!(xlabel="Time", ylabel="Amplitude", title=dataset)
				plot!(ylims=(-3, 3))
				savefig("./pipeline/plots/"*dataset*"/"*dataset*"_signal_"*segment*".png")

				# read entropy curves
				plot(size=(1024, 768))
				for r in r_list
					rcmse_curve = read(processed_file, "complexity/"*segment*"/RCMSE/"*string(r)*"/data")
					rcmse_scales = read(processed_file, "complexity/"*segment*"/RCMSE/"*string(r)*"/meta/scales")
					plot!(rcmse_scales, rcmse_curve, label="RCMSE_"*string(r))
				end
				plot!(xlabel="Scale", ylabel="Entropy", title=dataset)
				#y label from 0 to 3
				plot!(ylims=(0, 3))
				savefig("./pipeline/plots/"*dataset*"/"*dataset*"_rcmse_"*segment*".png")
			end
		end
	end
end


#### PLOT COMPLEXITY DISTRIBUTION ####

# linear regression slope dictionary
lr_values = Dict()
for group in g_list
	lr_values[group] = Dict()
	for r in r_list
		lr_values[group][r] = Dict()
		for sex in ["M", "F"]
			lr_values[group][r][sex] = Float64[]
		end
	end
end

# read linear regression values
for group in g_list
	for dataset in grouped_datasets[group]
		h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "r") do processed_file
			for r in r_list
				b = read(processed_file, "complexity/complete/RCMSE/"*r*"/meta/linear_fit/b")
				if dataset in sex["M"]
					push!(lr_values[group][r]["M"], b)
				elseif dataset in sex["F"]
					push!(lr_values[group][r]["F"], b)
				end
			end
		end
	end
end

println("Plotting complexity distribution...")
for r in r_list
	p = plot(size=(1024, 768))
	m_data = [lr_values[group][r]["M"] for group in g_list]
	f_data = [lr_values[group][r]["F"] for group in g_list]
	a_data = [[a;b] for (a, b) in zip(m_data, f_data)]

	m_labels = [l for (m, l) in zip(m_data, group_labels_plot) if length(m) > 0]
	m_labels = reshape(m_labels, 1, length(m_labels))
	m_data = [m for m in m_data if length(m) > 0]

	f_labels = [l for (f, l) in zip(f_data, group_labels_plot) if length(f) > 0]
	f_labels = reshape(f_labels, 1, length(f_labels))
	f_data = [f for f in f_data if length(f) > 0]

	violin!(m_labels, m_data, label="Male", line = 0, fill = (0.2, :green), side = :left)
	dotplot!(m_labels, m_data, label=false, line = 0, marker = (0.8, :green), side = :left)
	violin!(f_labels, f_data, label="Female", line = 0, fill = (0.2, :blue), side = :right)
	dotplot!(f_labels, f_data, label=false, line = 0, marker = (0.8, :blue), side = :right)
	boxplot!(group_labels_plot, a_data, label=false, line = (0.4, :black), fill = (0.2, :orange))
	plot!(xlabel="Group", ylabel="LRS", title="Linear regression slope distribution (RCMSE, 2, "*r*")")

	m_flag = false
	f_flag = false
	for x in p.series_list
		if x[:label] == "Male"
			if m_flag
				x[:label] = ""
			else
				m_flag = true
			end
		end
		if x[:label] == "Female"
			if f_flag
				x[:label] = ""
			else
				f_flag = true
			end
		end
	end

	savefig("./pipeline/plots/LRS_distribution_"*r*".png")
end

#plot complexity curves by group

for segment in s_list
	for r in r_list
		for group in g_list
			plot(size=(1024, 768))

			avg_curve = zeros(45)
			for dataset in grouped_datasets[group]
				h5open("./pipeline/processed_data/"*dataset*"/"*dataset*"_processed.h5", "r") do processed_file
					rcmse_curve = read(processed_file, "complexity/"*segment*"/RCMSE/"*string(r)*"/data")
					rcmse_scales = read(processed_file, "complexity/"*segment*"/RCMSE/"*string(r)*"/meta/scales")

					plot!(rcmse_scales, rcmse_curve, label=dataset)
					avg_curve += rcmse_curve
				end
			end
			avg_curve /= length(grouped_datasets[group])
			plot!(xlabel="Scale", ylabel="Entropy", title=group_labels[group])
			plot!(ylims=(0, 3))
			#plot average curve dotted line
			plot!([i for i in 1:45], avg_curve, label="Average", linestyle=:dash, line = (2, :black))
			savefig("./pipeline/plots/"*group_labels[group]*"_rcmse_"*string(r)*"_"*segment*".png")
		end
	end
end
