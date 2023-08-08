"""
This script computes the mean escape time as a function of the noise intensity for a range of
    different values of α.
"""


##################  Dependencies   ##################
using Makie, CairoMakie, Distributions, Random, ColorSchemes

################     Parameters      ################

nTransitions = 50
npoints = 5
npoints2 = 13
Δ = 2; 
ω = 2;
α = vcat([2;1.99;1.95], range(1.9,1,10))
λ = zeros(npoints,npoints2)
tippingType = zeros(npoints,1)
exponent = zeros(npoints,1)
initx = ω + rand()/10;

σ1 = range(2,0.6,npoints)
σ2 = range(2,0.1,npoints)

################# Main loop #####################

for i in 1:npoints2
    if i==1
        for j in 1:npoints
            transitions = getTransitions(nTransitions,initx,α[i],σ1[j],Δ,ω)
            λ[j,i], tippingType[j] = mean(transitions, dims = 1)
            print(100*((i-1)*npoints+j)/(npoints*npoints2))
            println("%")
        end
    elseif i > 1
        for j in 1:npoints
            transitions = getTransitions(nTransitions,initx,α[i],σ2[j],Δ,ω)
            λ[j,i], tippingType[j] = mean(transitions, dims = 1)
            print(100*((i-1)*npoints+j)/(npoints*npoints2))
            println("%")
        end
    end  
end

#################### Plotting #####################

mycolors = cgrad(:bluesreds,14; categorical = true)


kramersLaw =  ω^2*π/(Δ*sqrt(8)) * exp.(2Δ ./σ1.^2)
fig2 = lines(log.(σ1./Δ^(1/α[1])),log.(Δ*kramersLaw./ω.^α[1]), axis = (; limits = (-3,1.5,1,12), width = 500, height = 380, xlabel = "log (σ/Δ^1/α)", ylabel="log (Δ*λ/ω^α)"), label="Kramer's law",color = mycolors[1])

for i in 1:npoints2
    if i == 1
        lines!(log.(σ1./Δ^(1/α[i])),log.(Δ*λ[:,i]./ω.^α[i]), label = "α = $(α[i])",color = mycolors[i+1])
    else 
        lines!(log.(σ2./Δ^(1/α[i])),log.(Δ*λ[:,i]./ω.^α[i]), label = "α = $(α[i])",color = mycolors[i+1])
    end
end

axislegend()
fig2

save("trial alphaprofile.png",current_figure())
