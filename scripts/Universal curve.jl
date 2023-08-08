"""
This script computes the mean escape time as a function of the noise intensity for a range of
    different values of ω and Δ, and plots them the axis that define the universal curve.
"""


using Makie, CairoMakie, Distributions, Random, DrWatson, Tables, CSV, DataFrames, SpecialFunctions

################     Parameters      ################
nTransitions = 500;
npoints2 = 3;
npoints = 25;
Δ = range(1,3,npoints2);
ω = range(1,3,npoints2);
α = 1.95;
σ = range(3,0.2,npoints);

λ = zeros(npoints,npoints2,npoints2)

######################### Main loop #######################

for k in range(1,npoints2)
    initx = ω[k] + (rand() - 0.5)/10;
    for i in range(1,npoints2)
        for j in range(1,npoints)
            transitions = getTransitions(nTransitions,initx,α,σ[j],Δ[i],ω[k])
            λ[j,i,k] = mean(transitions[:,1])
            print(100*((k-1)*npoints2*npoints+(i-1)*npoints + j)/(npoints*npoints2^2))
            println("%")
        end
    end
end

###################### Plotting #######################
levyrate = (1.95 ./(√(2)*σ)).^α *2*gamma(1-α)*cos(π*α/2)
Kramerrate = π/√(8) * exp.(2 ./(σ .^2))

fig = lines(log.(σ./Δ[1]^(1/α)),log.(Δ[1]*λ[:,1,1]/ω[1]^α), axis = (; limits = (-2,1,0,10), width = 500, height = 380, xlabel="log (σ/Δ^1/α)", ylabel="log (Δ*λ/ω^α)"))#,label = "Δ = $(Δ[1]), ω = $(ω[1])");
lines!(log.(σ./Δ[2]^(1/α)),log.(Δ[2]*λ[:,2,1]/ω[1]^α))#,label = "Δ = $(Δ[2]), ω = $(ω[1])");
lines!(log.(σ./Δ[3]^(1/α)),log.(Δ[3]*λ[:,3,1]/ω[1]^α))#,label = "Δ = $(Δ[3]), ω = $(ω[1])");
for j in range(2,npoints2)
    for i in range(1,npoints2)
        lines!(log.(σ./Δ[i]^(1/α)),log.(Δ[i]*λ[:,i,j]/ω[j]^α))#,label = "Δ = $(Δ[i]), ω = $(ω[j])")
    end
end
lines!(log.(σ[1:23]./Δ[1]^(1/2)),log.(Δ[1]*Kramerrate[1:23]/ω[1]^2),label = "Kramer's time")
lines!(log.(σ./Δ[1]^(1/α)),log.(Δ[1]*levyrate/ω[1]^α),label = "Asymptotic escape time")
axislegend()
current_figure()

save("trial UniversalCurve.png", current_figure())

