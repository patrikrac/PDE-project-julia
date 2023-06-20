 #=
Simple timing functions for the EDP project
 =#
using edp
using LinearSolve

function f_rhs(x, b=(1,1), p=3, q=3, r=2)
    return exp((b[1]*x[1]+b[2]*x[2])/2.)*sin(p*r*pi*x[1])*sin(q*r*pi*x[2])
end

function alpha(b,p=3,q=3,r=2, c=1.)
    return 1.0/((b[1]^2 / 4.) + (b[2]^2 / 4.) + (pi^2 * r^2 * (p^2 + q^2)) + c)
end


function timing()
    h_base = 0.1
    h_values = [h_base * 2.0^(-i) for i = 0:6]

    for h ∈ h_values
        N = Int(1/h)
        println("h = $h (N = $N)")
        vtxL, eltL = GenerateLShapedMesh(N, Int(N/2))
        println("Number of vertices:\t$(length(vtxL))")
        t = @elapsed A, F = Assemble(vtxL, eltL, x -> f_rhs(x))
        println("Assemble:\t\t$(t)s")
        p = LinearProblem(A,F)
        t = @elapsed Uh = solve(p, KrylovJL_GMRES())
        println("Solve:\t\t\t$(t)s")
        println()
    end
end

timing()