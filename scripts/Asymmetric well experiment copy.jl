"""
This script computes the mean escape time of the potential wells as a function of the noise intensity for
    the asymmetric potential
"""



################### Dependencies ######################

using Makie, CairoMakie, Distributions, Random, DrWatson, Tables, CSV, SpecialFunctions, DataFrames

################     Parameters      ################

nTransitions = 500;
npoints = 50;

α = 1.9
σ = range(0.3,1,npoints)

λᵣ = zeros(npoints)
λₗ = zeros(npoints)
tippingType = zeros(npoints,1)

initx = 4 + rand()/30


##################### Main loop ################################

λright = zeros(npoints)
λleft = zeros(npoints)
for j in range(1,npoints)
    leftTransitions, rightTransitions = getTransitionsAsymmetric(nTransitions,initx,α,σ[j])
    λright[j], tippingType[j] = mean(rightTransitions, dims = 1) 
    λleft[j], tippingType[j] = mean(leftTransitions, dims = 1) 
    print(100*j/npoints)
    println("%")
end

################# Plotting ######################################

leftwellrate = (3.2 ./(√(2)*σ)).^α *2*gamma(1-α)*cos(π*α/2)

rightwellrate2 = (5.6./(√(2)*σ)).^α *2*gamma(1-α)*cos(π*α/2)
experimentalrate =  exp.(exp.(1.216) * σ.^(-0.759))

fig2 = lines(log.(σ),log.(λright), axis = (; title = "Transition time for α = $α", xlabel = "log σ", ylabel="log λ"), label="Right well")
lines!(log.(σ),log.(λleft), label = "Left well")
lines!(log.(σ),log.(rightwellrate2))
lines!(log.(σ),log.(leftwellrate))
lines!(log.(σ), log.(experimentalrate))
axislegend()
fig2

save("Trial Heterogeneous response.png",current_figure())


############ Save data ###########

result = DataFrame(Tables.table(hcat(σ,λright,λleft)))
rename!(result,:Column1 => :σ, :Column2 => :RightWell, :Column3=> :LeftWell)
params = @strdict α
CSV.write(savename("Trial Asymetric potential data ",params,"csv"),result)