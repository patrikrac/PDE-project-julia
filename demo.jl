
using edp
using LinearSolve
using SparseArrays

function f_rhs(x, b, p=3, q=3, r=2)
    return exp((b[1]*x[1]+b[2]*x[2])/2.)*sin(p*r*pi*x[1])*sin(q*r*pi*x[2])
end

function alpha(b,p,q,r, c=1.)
    return 1.0/((b[1]^2 / 4.) + (b[2]^2 / 4.) + (pi^2 * r^2 * (p^2 + q^2)) + c)
end

p = 3
q = 3
r = 2
b = (1, 1)
c = 1.

vtxL, eltL = GenerateLShapedMesh(1000, 500)
println("Generated Mesh")

α = alpha(b,p,q,r,c)

U = [α * f_rhs(x,b,p,q,r) for x in vtxL]

PlotApproximation(vtxL, eltL, U)

bndL = Boundary(eltL)
bndL_set = Set(v for e in bndL for v in e)

println("Generated Boundary")

M = MassMatrix(vtxL, eltL, bndL_set)
K = StiffnessMatrix(vtxL, eltL, bndL_set)
C = CentralMatrix(vtxL, eltL, bndL_set, b)

F = RightHandSide(vtxL, eltL, x -> f_rhs(x, b), bndL_set)



A = sparse(c.*M + K + C)
prob = LinearProblem(A,F)

Uh = solve(prob, KrylovJL_GMRES())
PlotApproximation(vtxL, eltL, Uh)

println("The error is $(L2(Uh-U, M) / L2(U,M))")

#=
M = globalMatrix(vtxL, eltL, Mloc)
U_test = ones(length(vtxL))

println(dot(U_test, M*U_test))

K = globalMatrix(vtxL, eltL, Kloc)

println(norm(K*U_test))
=#


vtx, elt = GenerateRectangleMesh(1,1,2,2)

display(vtx)
display(elt)

PlotMesh(vtx, elt)

bnd = Boundary(elt)
bnd_set = Set(v for e in bnd for v in e)

#=
using LinearSolve
using SparseArrays

A = sparse(rand(10, 10))
B = rand(10)
X = zeros(10)

prob = LinearProblem(A,B)

sol = solve(prob, KrylovJL_GMRES())

=#