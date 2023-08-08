"""
    getTransitions(nTransitions,initialCondition,α,σ,Δ,ω) --> transitions

Method to compute the transition times in a Lipschitz double well. This method takes the number the transitions 
wanted nTransitions, the initial condition of the process initialCondition, the stability parameter α, the noise 
amplitude σ, the width of the well ω and its height Δ; and returns an nTransitions × 2 array with the times between
one transition and the next, and a logic (1 or 0, not true or false) value that takes the value 1 if there was a
jump of size ω in the respective transition. This is ment to be a proxy for the transitions that occur through only
one large jump and those that cross the barrier in several steps. 
"""
function getTransitions(nTransitions,initialCondition,α,σ,Δ,ω)
    events = 0;
    time = 0.0;
    cross = 0.0;
    lastTime = 0.0;
    upperThr = ω - ω/10
    largestNoise = 0.0;
    trajectory = initialCondition;
    transitions = zeros(nTransitions,2)
    positive = initialCondition > 0;
    while events < nTransitions
        (newPoint, δt, largestNoise) = EMadaptivestep(trajectory,α,σ,Δ,ω,largestNoise)
        time += δt
        if newPoint == NaN
            error("Numerical Overflow D:")
        end
        if abs(newPoint) > upperThr
            if newPoint > upperThr
                if positive
                    cross = time;
                    largestNoise = 0.0;
                else
                    if events == 0
                        lastTime = time
                        cross = time
                        events += 1
                        transitions[events,1] = (time + cross)/2;
                        positive = true
                        if largestNoise > ω
                            transitions[events,2] = 1;
                        end
                        largestNoise = 0.0;
                    else
                        events += 1;
                        eventTime = (time + cross)/2;
                        transitions[events,1] = eventTime - lastTime;
                        lastTime = eventTime;
                        cross = time
                        positive = true;
                        if largestNoise > ω
                            transitions[events,2] = 1;
                        end
                        largestNoise = 0.0;
                    end                 
                end
            else
                if positive
                    if events == 0
                        lastTime = time
                        cross = time
                        events += 1
                        transitions[events,1] = (time + cross)/2;
                        positive = false
                        if largestNoise > ω
                            transitions[events,2] = 1;
                        end
                        largestNoise = 0.0;
                    else
                        events += 1;
                        eventTime = (time + cross)/2;
                        transitions[events,1] = eventTime - lastTime;
                        lastTime = eventTime;
                        cross = time
                        positive = false;
                        if largestNoise > ω
                            transitions[events,2] = 1;
                        end
                        largestNoise = 0.0;
                    end
                else
                    cross = time;  
                    largestNoise = 0.0;               
                end
            end
        end                
        trajectory = newPoint;
    end
    return transitions
end



"""
    getTransitionsNonLipschitz(nTransitions,initialCondition,α,σ,Δ,ω)

    Method to compute the transition times in a quartic double well. This method takes the number the transitions 
    wanted nTransitions, the initial condition of the process initialCondition, the stability parameter α, the noise 
    amplitude σ, the width of the well ω and its height Δ; and returns an nTransitions × 2 array with the times between
    one transition and the next, and a logic (1 or 0, not true or false) value that takes the value 1 if there was a
    jump of size ω in the respective transition. This is ment to be a proxy for the transitions that occur through only
    one large jump and those that cross the barrier in several steps. 
"""
function getTransitionsNonLipschitz(nTransitions,initialCondition,α,σ,Δ,ω)
    events = 0;
    time = 0.0;
    cross = 0.0;
    lastTime = 0.0;
    upperThr = ω - ω/10
    largestNoise = 0.0;
    trajectory = initialCondition;
    transitions = zeros(nTransitions,2)
    positive = initialCondition > 0;
    while events < nTransitions
        (newPoint, δt, largestNoise) = EMadaptivestepNonLipschitz(trajectory,α,σ,Δ,ω,largestNoise)
        time += δt
        if newPoint == NaN
            error("Numerical Overflow D:")
        end
        if abs(newPoint) > upperThr
            if newPoint > upperThr
                if positive
                    cross = time;
                    largestNoise = 0.0;
                else
                    if events == 0
                        lastTime = time
                        cross = time
                        events += 1
                        transitions[events,1] = (time + cross)/2;
                        positive = true
                        if largestNoise > ω
                            transitions[events,2] = 1;
                        end
                        largestNoise = 0.0;
                    else
                        events += 1;
                        eventTime = (time + cross)/2;
                        transitions[events,1] = eventTime - lastTime;
                        lastTime = eventTime;
                        cross = time
                        positive = true;
                        if largestNoise > ω
                            transitions[events,2] = 1;
                        end
                        largestNoise = 0.0;
                    end                 
                end
            else
                if positive
                    if events == 0
                        lastTime = time
                        cross = time
                        events += 1
                        transitions[events,1] = (time + cross)/2;
                        positive = false
                        if largestNoise > ω
                            transitions[events,2] = 1;
                        end
                        largestNoise = 0.0;
                    else
                        events += 1;
                        eventTime = (time + cross)/2;
                        transitions[events,1] = eventTime - lastTime;
                        lastTime = eventTime;
                        cross = time
                        positive = false;
                        if largestNoise > ω
                            transitions[events,2] = 1;
                        end
                        largestNoise = 0.0;
                    end
                else
                    cross = time;  
                    largestNoise = 0.0;               
                end
            end
        end                
        trajectory = newPoint;
    end
    return transitions
end




"""
    getTransitionsAsymmetric(nTransitions,initialCondition,α,σ)

    Method to compute the transition times in the asymmetric double well. This method takes the number the transitions 
    wanted nTransitions, the initial condition of the process initialCondition, the stability parameter α, the noise 
    amplitude σ, the width of the well ω and its height Δ; and returns an nTransitions × 2 array with the times between
    one transition and the next, and a logic (1 or 0, not true or false) value that takes the value 1 if there was a
    jump of size ω in the respective transition. This is ment to be a proxy for the transitions that occur through only
    one large jump and those that cross the barrier in several steps. 
"""
function getTransitionsAsymmetric(nTransitions,initialCondition,α,σ)
    eventsleft = 0;
    eventsright = 0;
    time = 0.0;
    cross = 0.0;
    lastTime = 0.0;
    totalTransitions = 2 * nTransitions
    upperThr1 = -0.9
    upperThr2 = 2.7
    largestNoise = 0.0;
    trajectory = initialCondition;
    leftTransitions = zeros(nTransitions,2)
    rightTransitions = zeros(nTransitions,2)
    positive = initialCondition > 0;
    while eventsleft + eventsright < totalTransitions
        (newPoint, δt, largestNoise) = EMadaptivestepAsymmetric(trajectory,α,σ,largestNoise)
        time += δt
        if newPoint == NaN
            error("Numerical Overflow D:")
        end
        if newPoint > upperThr2
            if positive
                cross = time;
                largestNoise = 0.0;
            else
                if eventsleft == 0
                    lastTime = time
                    cross = time
                    eventsleft += 1
                    leftTransitions[eventsleft,1] = (time + cross)/2;
                    positive = true
                    if largestNoise > 1
                        leftTransitions[eventsleft,2] = 1;
                    end
                    largestNoise = 0.0;
                else
                    eventsleft += 1;
                    eventTime = (time + cross)/2;
                    leftTransitions[eventsleft,1] = eventTime - lastTime;
                    lastTime = eventTime;
                    cross = time
                    positive = true;
                    if largestNoise > 1
                        leftTransitions[eventsleft,2] = 1;
                    end
                    largestNoise = 0.0;
                end                 
            end
        elseif newPoint < upperThr1
            if positive
                if eventsright == 0
                    lastTime = time
                    cross = time
                    eventsright += 1
                    rightTransitions[eventsright,1] = (time + cross)/2;
                    positive = false
                    if largestNoise > 2
                        rightTransitions[eventsright,2] = 1;
                    end
                    largestNoise = 0.0;
                else
                    eventsright += 1;
                    eventTime = (time + cross)/2;
                    rightTransitions[eventsright,1] = eventTime - lastTime;
                    lastTime = eventTime;
                    cross = time
                    positive = false;
                    if largestNoise > 2
                        rightTransitions[eventsright,2] = 1;
                    end
                    largestNoise = 0.0;
                end
            else
                cross = time;  
                largestNoise = 0.0;               
            end
        end 
        trajectory = newPoint;
    end
    return [leftTransitions,rightTransitions]
end



"""
    getTransitionsAsymmetricLeft(nTransitions,α,σ)

    Method to compute the transition times from the left well of the asymmetric potential. For debugging and comparing 
    against. This method takes the number the transitions wanted nTransitions, the initial condition of the process 
    initialCondition, the stability parameter α, the noise amplitude σ, the width of the well ω and its height Δ; and
    returns an nTransitions × 2 array with the times between one transition and the next, and a logic (1 or 0, not true
    or false) value that takes the value 1 if there was a jump of size ω in the respective transition. This is ment to
    be a proxy for the transitions that occur through only one large jump and those that cross the barrier in several steps. 
"""
function getTransitionsAsymmetricLeft(nTransitions,α,σ)
    events = 0;
    time = 0.0;
    upperThr = 2.7
    largestNoise = 0.0;
    initialCondition = rand()/20 - 1;
    trajectory = initialCondition;
    leftTransitions = zeros(nTransitions,2)
    while events < nTransitions
        (newPoint, δt, largestNoise) = EMadaptivestepAsymmetric(trajectory,α,σ,largestNoise)
        time += δt
        if newPoint == NaN
            error("Numerical Overflow D:")
        end
        if newPoint > upperThr
            events += 1
            leftTransitions[events,1] = time;
            if largestNoise > 1
                leftTransitions[events,2] = 1;
            end
            largestNoise = 0.0;
            time = 0.0;
            trajectory = initialCondition
        else
            trajectory = newPoint;
        end
    end
    return leftTransitions
end


"""
    getTransitionsAsymmetricLeft(nTransitions,α,σ)

    Method to compute the transition times from the right well of the asymmetric potential. For debugging and comparing 
    against. This method takes the number the transitions wanted nTransitions, the initial condition of the process 
    initialCondition, the stability parameter α, the noise amplitude σ, the width of the well ω and its height Δ; and
    returns an nTransitions × 2 array with the times between one transition and the next, and a logic (1 or 0, not true
    or false) value that takes the value 1 if there was a jump of size ω in the respective transition. This is ment to
    be a proxy for the transitions that occur through only one large jump and those that cross the barrier in several steps. 
"""
function getTransitionsAsymmetricRight(nTransitions,α,σ)
    events = 0;
    time = 0.0;
    upperThr = -0.9
    largestNoise = 0.0;
    initialCondition = rand()/20 + 3;
    trajectory = initialCondition;
    rightTransitions = zeros(nTransitions,2)
    while events < nTransitions
        (newPoint, δt, largestNoise) = EMadaptivestepAsymmetric(trajectory,α,σ,largestNoise)
        time += δt
        if newPoint == NaN
            error("Numerical Overflow D:")
        end
        if newPoint < upperThr
            events += 1
            rightTransitions[events,1] = time;
            if largestNoise > 1
                rightTransitions[events,2] = 1;
            end
            largestNoise = 0.0;
            time = 0.0;
            trajectory = initialCondition
        else
            trajectory = newPoint;
        end
    end
    return rightTransitions
end
