
"""
This file contains the methods that return an integrated time step of the Langevin equation with the respective
potential well. If you want to use another potential, write a method that returns the slope in potentialDerivative.jl
and modify the next template changing the three instances of your_function(param), and the name of the function 


function YourEMadaptivestep(lastValue,α,σ,Δ,ω,largestNoiseSoFar)
    minTimestep = 0.0001;
    maxTimestep = 0.01;
    errorTolerance = 0.0001;
    δt = max(minTimestep,min(maxTimestep,1/abs(your_function(param)) * (errorTolerance + sqrt(errorTolerance^2 + 4 * errorTolerance))/2));
    v =π*( -0.5 + rand());
    w =-log(rand());
    η =(sin(α*v)/cos(v)^(1/α))*(cos((1-α)*v)/w)^((1-α)/α);
    if δt > minTimestep

        noise = σ/√2 * δt^(1/α) * η

        if abs(noise) > largestNoiseSoFar
            largestNoiseSoFar = abs(noise);
            #println(largestNoise)
        end
        step = lastValue - δt*(your_function(param)) + noise;
        
    else
        pot = your_function(param);
        noise = σ/√2 * δt^(1/α) * η
        if abs(noise) > largestNoiseSoFar
            largestNoiseSoFar = abs(noise);
        end
        step = lastValue - δt*(pot/(1 + δt * abs(pot))) + noise;
        
    end
    return step, δt, largestNoiseSoFar
end

"""


"""
    EMadaptivestep(lastValue,α,σ,Δ,ω,largestNoiseSoFar) --> step, δt, largestNoiseSoFar

Integrates one timestep using the adaptive tamed Euler-Maruyama scheme from the paper
ADAPTIVE TIMESTEPPING STRATEGIES FOR NONLINEAR STOCHASTIC SYSTEMS∗ CONALL KELLY  AND GABRIEL J. LORD.
The inputs are the last state lastValue, the stability parameter α, the noise amplitude σ, the width of the well ω, 
the height Δ and the largest noise observed so far in the process largestNoiseSoFar.
The method outputs the next state step, the length of the time step δt and the the largest noise observed so far 
in the process largestNoiseSoFar.
"""
function EMadaptivestep(lastValue,α,σ,Δ,ω,largestNoiseSoFar)
    minTimestep = 0.0001;
    maxTimestep = 0.01;
    errorTolerance = 0.0001;
    δt = max(minTimestep,min(maxTimestep,1/abs(potentialDerivative(lastValue,Δ,ω)) * (errorTolerance + sqrt(errorTolerance^2 + 4 * errorTolerance))/2));
    v =π*( -0.5 + rand());
    w =-log(rand());
    η =(sin(α*v)/cos(v)^(1/α))*(cos((1-α)*v)/w)^((1-α)/α);
    if δt > minTimestep

        noise = σ/√2 * δt^(1/α) * η

        if abs(noise) > largestNoiseSoFar
            largestNoiseSoFar = abs(noise);
            #println(largestNoise)
        end
        step = lastValue - δt*(potentialDerivative(lastValue,Δ,ω)) + noise;
        
    else
        pot = potentialDerivative(lastValue,Δ,ω);
        noise = σ/√2 * δt^(1/α) * η
        if abs(noise) > largestNoiseSoFar
            largestNoiseSoFar = abs(noise);
        end
        step = lastValue - δt*(pot/(1 + δt * abs(pot))) + noise;
        
    end
    return step, δt, largestNoiseSoFar
end




"""
    EMadaptivestepNonLipschitz(lastValue,α,σ,Δ,ω,largestNoiseSoFar) --> step, δt, largestNoiseSoFar

Integrates one timestep using the adaptive tamed Euler-Maruyama scheme from the paper
ADAPTIVE TIMESTEPPING STRATEGIES FOR NONLINEAR STOCHASTIC SYSTEMS∗ CONALL KELLY  AND GABRIEL J. LORD.
The inputs are the last state lastValue, the stability parameter α, the noise amplitude σ, the width of the well ω, 
the height Δ and the largest noise observed so far in the process largestNoiseSoFar.
The method outputs the next state step, the length of the time step δt and the the largest noise observed so far 
in the process largestNoiseSoFar. Note that this function employs a different potential, see potentialDerivativeNonLipschitz.
"""
function EMadaptivestepNonLipschitz(lastValue,α,σ,Δ,ω,largestNoiseSoFar)
    minTimestep = 0.0001;
    maxTimestep = 0.01;
    errorTolerance = 0.0001;
    δt = max(minTimestep,min(maxTimestep,1/abs(potentialDerivativeNonLipschitz(lastValue,Δ,ω)) * (errorTolerance + sqrt(errorTolerance^2 + 4 * errorTolerance))/2));
    v =π*( -0.5 + rand());
    w =-log(rand());
    η =(sin(α*v)/cos(v)^(1/α))*(cos((1-α)*v)/w)^((1-α)/α);
    if δt > minTimestep
        noise = σ/√2 * δt^(1/α) * η
        if abs(noise) > largestNoiseSoFar
            largestNoiseSoFar = abs(noise);
            #println(largestNoise)
        end
        step = lastValue - δt*(potentialDerivativeNonLipschitz(lastValue,Δ,ω)) + noise;
        
    else
        pot = potentialDerivativeNonLipschitz(lastValue,Δ,ω);
        noise = σ/√2 * δt^(1/α) * η
        if abs(noise) > largestNoiseSoFar
            largestNoiseSoFar = abs(noise);
        end
        step = lastValue - δt*(pot/(1 + δt * abs(pot))) + noise;
        
    end
    return step, δt, largestNoiseSoFar
end




"""
    EMadaptivestepAsymmetric(lastValue,α,σ,largestNoiseSoFar) --> step, δt, largestNoiseSoFar

Integrates one timestep using the adaptive tamed Euler-Maruyama scheme from the paper
ADAPTIVE TIMESTEPPING STRATEGIES FOR NONLINEAR STOCHASTIC SYSTEMS∗ CONALL KELLY  AND GABRIEL J. LORD.
The inputs are the last state lastValue, the stability parameter α, the noise amplitude σ, the width of the well ω, 
the height Δ and the largest noise observed so far in the process largestNoiseSoFar.
The method outputs the next state step, the length of the time step δt and the the largest noise observed so far 
in the process largestNoiseSoFar. Note that this function employs a different potential, see potentialDerivativeAsymmetric.
"""
function EMadaptivestepAsymmetric(lastValue,α,σ,largestNoiseSoFar)
    minTimestep = 0.0001;
    maxTimestep = 0.01;
    errorTolerance = 0.0001;
    pot = potentialDerivativeAsymmetric(lastValue);
    δt = max(minTimestep,min(maxTimestep,1/abs(pot) * (errorTolerance + sqrt(errorTolerance^2 + 4 * errorTolerance))/2));
    v =π*( -0.5 + rand());
    w =-log(rand());
    η =(sin(α*v)/cos(v)^(1/α))*(cos((1-α)*v)/w)^((1-α)/α);
    if δt > minTimestep
        noise = σ/√2 * δt^(1/α) * η
        if abs(noise) > largestNoiseSoFar
            largestNoiseSoFar = abs(noise);
            #println(largestNoise)
        end
        step = lastValue - δt*(pot) + noise;
        
    else
        noise = σ/√2 * δt^(1/α) * η
        if abs(noise) > largestNoiseSoFar
            largestNoiseSoFar = abs(noise);
        end
        step = lastValue - δt*(pot/(1 + δt * abs(pot))) + noise;
        
    end
    return step, δt, largestNoiseSoFar
end
