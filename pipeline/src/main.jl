# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

include("processing.jl")

filter = Dict(
	"type" => "BPF",
	"cutoff_low" => 0.1, # Hz
	"cutoff_high" => 100 # Hz
)

groups = ["A", "B", "C", "D", "E", "F", "G", "H"]

grouped_datasets = Dict()

grouped_datasets["A"] = ["MR-0311", "MR-0309", "MR-0306", "MR-0303", "MR-0300-t1", "MR-0299-t2", "MR-0298-t1", "MR-0296-t2"]
grouped_datasets["B"] = ["MR-0294", "MR-0289", "MR-0288-t1", "MR-0284", "MR-0283-t1", "MR-0282", "MR-0276", "MR-0273", "MR-0270"]
grouped_datasets["C"] = ["MR-0313", "MR-0312", "MR-0310", "MR-0307-t2", "MR-0305", "MR-0304-t2", "MR-0302-t1", "MR-0301-t2", "MR-0297"]
grouped_datasets["D"] = ["MR-0293", "MR-0292-t1", "MR-0291", "MR-0290", "MR-0287", "MR-0285-t1", "MR-0280-t1", "MR-0278", "MR-0275", "MR-0274"]
grouped_datasets["E"] = ["MR-0593_nd4", "MR-0592_nd4", "MR-0591_nd4", "MR-0573_nd4", "MR-0483", "MR-0460", "MR-0456"]
grouped_datasets["F"] = ["MR-0599_nd4", "MR-0597_nd4", "MR-0596_nd4", "MR-0569_nd4", "MR-0554", "MR-0465-t2"]
grouped_datasets["G"] = ["MR-0586_nd4", "MR-0585_nd4", "MR-0583_nd4", "MR-0579_nd4", "MR-0577_nd4", "MR-0575_nd4"]
grouped_datasets["H"] = ["MR-0588_nd4", "MR-0587_nd4", "MR-0582_nd4", "MR-0568_nd4", "MR-0552", "MR-0548-t1", "MR-0530"]

datasets = []
for group in groups
	global datasets = [datasets; grouped_datasets[group]]
end

# parameters
snr_threshold = 10
resampling_rate = 250 # Hz
r = 0.2

for dataset in datasets
	println("Processing dataset: ", dataset)
	# if directory does not exist, create it
	if !isdir("./pipeline/processed_data/"*dataset*"/")
		mkpath("./pipeline/processed_data/"*dataset*"/")
	end

	# processing pipeline
	select_with_signal_to_noise_ratio(dataset)
	normalize_signals_and_average(dataset)
	filter_signal(dataset, filter)
	resample_signal(dataset, resampling_rate)

	compute_complexity_curve(dataset, "MSE", 2, r, [i for i in 1:45])
	compute_complexity_curve(dataset, "RCMSE", 2, r, [i for i in 1:45])
	compute_complexity_curve(dataset, "FMSE", 2, r, [i for i in 1:45])
	compute_complexity_curve(dataset, "FRCMSE", 2, r, [i for i in 1:45])
end