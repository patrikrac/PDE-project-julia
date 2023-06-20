module edp

include("geometry.jl")
include("utilities.jl")
include("assembly.jl")

export GenerateRectangleMesh, GenerateLShapedMesh

export PlotMesh, PlotApproximation, Boundary, L2

export MassMatrix, StiffnessMatrix, CentralMatrix, RightHandSide, Assemble

end # module edp
