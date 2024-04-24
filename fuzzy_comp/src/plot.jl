# -*- coding: utf-8 -*-
# @author: max (max@mdotmonar.ch)

using Plots

#open txt
f = open("./fuzzy_comp/processed_data/fmse_noise_julia.txt")
fmse_s = read(f, String)
close(f)
#parse txt
fmse_s = split(fmse_s, "\n")
fmse_julia = Float64[]
for n in fmse_s
    if n != ""
        push!(fmse_julia, parse(Float64, n))
    end
end
println(fmse_julia)

#open txt
f = open("./fuzzy_comp/processed_data/mse_noise_julia.txt")
mse_s = read(f, String)
close(f)
#parse txt
mse_s = split(mse_s, "\n")
mse_julia = Float64[]
for n in mse_s
    if n != ""
        push!(mse_julia, parse(Float64, n))
    end
end
println(mse_julia)

#open txt
f = open("./fuzzy_comp/processed_data/fmse_noise_python.txt")
fmse_s = read(f, String)
close(f)
#parse txt
fmse_s = split(fmse_s, "\n")
fmse_python = Float64[]
for n in fmse_s
    if n != ""
        push!(fmse_python, parse(Float64, n))
    end
end
println(fmse_python)

#open txt
f = open("./fuzzy_comp/processed_data/mse_noise_python.txt")
mse_s = read(f, String)
close(f)
#parse txt
mse_s = split(mse_s, "\n")
mse_python = Float64[]
for n in mse_s
    if n != ""
        push!(mse_python, parse(Float64, n))
    end
end
println(mse_python)

#plot
plot(size=(800, 600))
plot!(mse_julia, label="MSE Julia")
plot!(fmse_julia, label="FMSE Julia")
plot!(mse_python, label="MSE Python")
plot!(fmse_python, label="FMSE Python")
plot!(xlabel="Scale", ylabel="Entropy", title="Noise")
plot!(ylims=(0, 1))
savefig("./fuzzy_comp/plots/noise.png")