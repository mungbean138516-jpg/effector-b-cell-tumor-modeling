using Pkg

# Run from the repository root.
Pkg.activate(".")

Pkg.add("BifurcationKit")
Pkg.add("Accessors")
Pkg.add("ForwardDiff")
Pkg.instantiate()

println("Two-parameter continuation environment is ready.")
println("Now run:")
println("include(\"continue_saddle_node_in_b5_b6.jl\")")
