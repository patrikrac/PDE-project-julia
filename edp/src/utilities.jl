#=
Implemenation of some utility functions
 -> Plotting of the mesh
=#

using Meshes, MeshViz
using ColorSchemes
using LinearAlgebra

import GLMakie as Mke

function PlotMesh(vtx, elt)
    mesh = SimpleMesh(vtx, connect.(elt))
    viz(mesh, showfacets = true)
end

function PlotMesh(vtx, elt, val)
    mesh = SimpleMesh(vtx, connect.(elt))
    viz(mesh, showfacets = true, color = val, colormap = :jet)
end

function PlotApproximation(vtx, elt, val)
    mesh = SimpleMesh(vtx, connect.(elt))
    viz(mesh, showfacets = false, color = val, colormap = :jet)
end


#=
Function to compute the boundary of a mesh
=#
function Boundary(elt)
    b = Set{Tuple{Int,Int}}()
    bnd_set = Set{Tuple{Int,Int}}()

    for tri ∈ elt
        for i ∈ 1:3
            edge = (tri[i], tri[mod(i,3)+1])
            if !(reverse(edge) ∈ bnd_set)
                push!(b, edge)
                push!(bnd_set, edge)
            else
                pop!(b, reverse(edge))
            end
        end
    end

    return collect(b)
end


#=
Generally usefull functions
=#
function L2(v, M)
    return sqrt(dot(v, M*v))
end
