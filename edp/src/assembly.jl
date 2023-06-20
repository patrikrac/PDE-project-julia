#=
FIle implementing of the matrix associated to the linrar system of equations
=#
using LinearAlgebra
using SparseArrays

#=
# Local matrix assembly
=#

# Mass matrix
function Mloc(vtx, e)
    v1 = vtx[e[1]]
    v2 = vtx[e[2]]
    v3 = vtx[e[3]]

    τ = abs((v1[1] - v3[1]) * (v2[2] - v1[2]) - (v1[1] - v2[1]) * (v3[2] - v1[2])) / 2.

    return τ/12. * [2. 1. 1.; 1. 2. 1.; 1. 1. 2.]
end

# Stiffness matrix
function Kloc(vtx, e)
    v1 = vtx[e[1]]
    v2 = vtx[e[2]]
    v3 = vtx[e[3]]

    v = [collect(v1), collect(v2), collect(v3)]

    τ = abs((v1[1] - v3[1]) * (v2[2] - v1[2]) - (v1[1] - v2[1]) * (v3[2] - v1[2])) / 2.

    K = zeros(3,3)
    for i ∈ 1:3
        for j = 1:3
            K[i,j] = ⋅(v[mod(i, 3) + 1] - v[mod(i+1, 3) + 1], v[mod(j, 3) + 1] - v[mod(j+1, 3) + 1])
        end
    end
    
    return K / (4*τ)
end

#Central matrix
function Cloc(vtx, e, b = (1,1))
    v1 = vtx[e[1]]
    v2 = vtx[e[2]]
    v3 = vtx[e[3]]

    v = [collect(v1), collect(v2), collect(v3)]

    τ = abs((v1[1] - v3[1]) * (v2[2] - v1[2]) - (v1[1] - v2[1]) * (v3[2] - v1[2])) / 2.

    C = zeros(3,3)
    for k ∈ 1:3
        nk = ×([0;0;1], [v[mod(k, 3) + 1] - v[mod(k+1, 3) + 1];0])[1:end-1]
        C[:, k] .= ⋅(b, nk) /  ⋅(v[k] - v[mod(k,3)+1], nk)
    end

    return τ / 3. * C
end


# Right hand side
function Floc(vtx, e, f)
    v1 = vtx[e[1]]
    v2 = vtx[e[2]]
    v3 = vtx[e[3]]

    v = [collect(v1), collect(v2), collect(v3)]

    τ = abs((v1[1] - v3[1]) * (v2[2] - v1[2]) - (v1[1] - v2[1]) * (v3[2] - v1[2])) / 2.

    F = zeros(3)
    for i ∈ 1:3
        ni = ×([0;0;1], [v[mod(i, 3) + 1] - v[mod(i+1, 3) + 1];0])[1:end-1]
        g(x) = f(x) * ⋅(x-v[mod(i,3)+1], ni) /  ⋅(v[i] - v[mod(i,3)+1], ni) 
        F[i] = g((v[1] + v[2]) / 2.) + g((v[2] + v[3]) / 2.) + g((v[3] + v[1]) / 2.)
    end

    return τ/3. * F
end



#=
# Global matrix assembly
=#
function globalMatrix(vtx, elt, locMatrix)
    nbr_vtx = length(vtx)
    nbr_elt = length(elt)

    d = length(elt[1])

    V = zeros(d, d, nbr_elt)
    I = zeros(Int, d, d, nbr_elt)
    J = zeros(Int, d, d, nbr_elt)

    for i ∈ 1:nbr_elt
        V[:,:,i] .= locMatrix(vtx, elt[i])
        I[:,:,i] .= elt[i]
        J[:,:,i] .= elt[i]
    end

    return sparse(I[:], J[:], V[:])
end


function globalMatrix(vtx, elt, bnd_vtx, locMatrix)
    nbr_vtx = length(vtx)
    nbr_elt = length(elt)

    d = length(elt[1])

    V = zeros(d, d, nbr_elt)
    I = zeros(Int, d, d, nbr_elt)
    J = zeros(Int, d, d, nbr_elt)

    for i ∈ 1:nbr_elt
        loc = locMatrix(vtx, elt[i])
        for j ∈ 1:d
            for k ∈ 1:d
                I[j,k,i] = elt[i][j]
                J[j,k,i] = elt[i][k]
                if (elt[i][j] ∉ bnd_vtx) && (elt[i][k] ∉ bnd_vtx)
                    V[j,k,i] = loc[j,k]
                elseif elt[i][j] == elt[i][k]
                    V[j,k,i] = 1.
                end
            end
        end
    end

    return sparse(I[:], J[:], V[:])
end


function MassMatrix(vtx, elt)
    return globalMatrix(vtx, elt, Mloc)
end

function MassMatrix(vtx, elt, bnd_vtx)
    return globalMatrix(vtx, elt, bnd_vtx, Mloc)
end


function StiffnessMatrix(vtx, elt)
    return globalMatrix(vtx, elt, Kloc)
end

function StiffnessMatrix(vtx, elt, bnd_vtx)
    return globalMatrix(vtx, elt, bnd_vtx, Kloc)
end

function CentralMatrix(vtx, elt, b)
    return globalMatrix(vtx, elt, Cloc)
end

function CentralMatrix(vtx, elt, bnd_vtx, b)
    return globalMatrix(vtx, elt, bnd_vtx, Cloc)
end


function RightHandSide(vtx, elt, f, bnd_vtx = Set())
    nbr_vtx = length(vtx)
    
    F = zeros(nbr_vtx)

    for e ∈ elt
        F_loc = Floc(vtx, e, f)
        for i ∈ 1:3
            if e[i] ∉ bnd_vtx
                F[e[i]] += F_loc[i]
            end
        end
    end

    return F
end


function Assemble(vtx, elt, f, b=(1,1), c=1.)
    bnd_vtx = Set(v for e in Boundary(elt) for v in e)

    M = MassMatrix(vtx, elt, bnd_vtx)
    K = StiffnessMatrix(vtx, elt, bnd_vtx)
    C = CentralMatrix(vtx, elt, bnd_vtx, b)

    F = RightHandSide(vtx, elt, f, bnd_vtx)

    A = sparse(c.*M + K + C)

    return A, F
end