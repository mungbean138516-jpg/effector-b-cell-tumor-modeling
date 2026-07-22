using Pkg

# Run this once from the repository root.
Pkg.activate(".")
Pkg.add("BifurcationKit")
Pkg.add("ForwardDiff")
Pkg.instantiate()

println("Bifurcation environment is ready.")
println("Now run:")
println("include(\"b6_formal_palc_continuation.jl\")")
