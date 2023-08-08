using DrWatson
@quickactivate "Alpha-transitions"

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())

Have fun with your new project!

You can help us improve DrWatson by opening
issues on GitHub, submitting feature requests,
or even opening your own Pull Requests!
"""
)
include("src//potentialDerivative.jl");
include("src//EMadaptivestep.jl");
include("src//getTransitions.jl");

println("Source code compiled!")