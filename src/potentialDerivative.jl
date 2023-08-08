
"""
This file contains the methods that return a the slope of the potential well. Note that the potentials are hardcoded,
if you want to use another potenial, write your method and follow the instructrions in EMadaptivestep.jl 
to integrate steps with it.

"""


"""
    potentialDerivative(x, Δ, ω, changeThreshold, limitSlope) --> drift
Returns the value of the drift, the derivative of the double well potential. The potential is made Lipschitz 
    by a smooth continuation that turns the quartic polinomial into a straight line of slope ±limitSlope in ±∞.
    Note that the existenceCondition imposes a relationship between the variables for the potential to exist:
        0 ≤ (2*changeThreshold*(limitSlope + 4*Δ/ω^2*changeThreshold - 4*Δ/ω^4*changeThreshold^3)-1)^2 - 
        4*changeThreshold^2*(limitSlope + 4*Δ/ω^2*changeThreshold - 4*Δ/ω^4*changeThreshold^3)^2)
        
"""
function potentialDerivative(x, Δ, ω, changeThreshold, limitSlope)
    auxCte1 = -4*Δ*changeThreshold/ω^2 + 4*Δ*changeThreshold^3/ω^4;
    existenceCondition = (2*changeThreshold*(limitSlope - auxCte1)-1)^2 - 4*changeThreshold^2*(limitSlope - auxCte1)^2;
    if existenceCondition < 0
        error("Modify parameters! The quantity (2*changeThreshold*(limitSlope + 4*Δ/ω^2*changeThreshold - 4*Δ/ω^4*changeThreshold^3)-1)^2-4*changeThreshold^2*(limitSlope + 4*Δ/ω^2*changeThreshold - 4*Δ/ω^4*changeThreshold^3)^2) must be positive and evaluates to $existenceCondition")
    end
    auxCte2 = (-2 * changeThreshold * (limitSlope - auxCte1)+1+sqrt(existenceCondition))/(2*changeThreshold^2*(limitSlope - auxCte1)+1)
    if abs(x) < changeThreshold
        -4*Δ/ω^2*x + 4*Δ/ω^4*x^3
    elseif x < 0
        -limitSlope+auxCte2/(1-auxCte2*x)^2   
    else
        limitSlope-auxCte2/(auxCte2*x+1)^2
    end
end





"""
    potentialDerivative(x, Δ, ω) --> drift
    Returns the value of the drift, the derivative of the double well potential. The potential is made Lipschitz 
    so the that turns the quartic polinomial into a straight line of slope ±limitSlope in ±∞.
        
"""

function potentialDerivative(x,Δ,ω)
    spikyparameter = 200
    Δ*x*(ω^2 - x^2)/(ω^4*(spikyparameter^4 + (x^2 - ω^2)^2 / ω^4)^(3/4)*(spikyparameter - (spikyparameter^4 + 1)^(1 / 4)))
end





"""
    potentialDerivativeNonLipschitz(x, Δ, ω) --> drift
    Returns the value of the drift, the derivative of the double well potential. The potential is not Lipschitz,
    is a quartic polinomial.
        
"""

function potentialDerivativeNonLipschitz(x,Δ,ω)
    4*Δ*x*(x^2 - ω^2)/ω^4
end





"""
    potentialDerivativeNonDifferentiable(x, Δ, ω) --> drift

Returns the value of the drift, the derivative of the double well potential. The potential is not differentiable,
    is piecewise linear.
        
"""

function potentialDerivativeNonDifferentiable(x,Δ,ω)
    if x <= -ω
        10
    elseif x < 0
        -Δ/ω
    elseif x < ω
        Δ/ω
    else
        -10
    end
end




"""
    potentialDerivativeAsymmetric2(x) --> drift

Returns the derivative of an asymmetric potential given by a polinomial of degree 6, auxiliar function.
    
"""
function potentialDerivativeAsymmetric2(x)
    (15 * x^5 + 135/8 * x^4 - 9 * x^3 - 135/8 * x^2 - 6 * x)
end



"""
    potentialDerivativeAsymmetric(x) --> drift

Returns the derivative of an asymmetric potential given by a polinomial of degree 6, 
    evaluated in a non-linear transformation of the x axis.
"""
function potentialDerivativeAsymmetric(x)
    potentialDerivativeAsymmetric2(2(tanh(-x/2-1)+1)-x/4) * (-(1-tanh(-x/2-1)^2)-1/4) 
end