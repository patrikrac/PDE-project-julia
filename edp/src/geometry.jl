#=
Geometry module for Partial differential equations module
    -> Contains functions to create rectangular and L-shaped domains
=#

function GenerateRectangleMesh(Lx, Ly, Nx, Ny)
    
    # Generate the vertices of the mesh based on the given parameters
    vtx = [(xv,yv) for yv in LinRange(0,Ly,Ny+1) for xv in LinRange(0,Lx,Nx+1)]

    elt = []

    # Generate the elements of the mesh based on the given parameters
    for iy=0:Ny-1
        for ix=0:Nx-1
             push!(elt, (iy*(Nx+1)+ix+1, iy*(Nx+1)+ix+2, (iy+1)*(Nx+1)+ix+2), (iy*(Nx+1)+ix+1, (iy+1)*(Nx+1)+ix+2, (iy+1)*(Nx+1)+ix+1))
        end
    end

    return vtx, elt
end


function GenerateLShapedMesh(N, Nl)
    h = 1.0/N
    l = h*Nl

    vtxl, eltl = GenerateRectangleMesh(1, l, N, Nl)

    if l == 1
        return vtxl, eltl
    end

    vtxu, eltu = GenerateRectangleMesh(l, 1-l, Nl, Int((1-l)/h))

    vtxu = [(xv,yv+l) for (xv,yv) in vtxu]

    vtx = union(vtxl, vtxu)

    interface_dict = Dict(v=>i for (i,v) in enumerate(vtxl))

    offset = length(vtxl) - (Nl+1)

    elt = Vector{Tuple{Int,Int,Int}}(undef, length(eltl) + length(eltu))

    for (i,e)∈enumerate(eltl)
        elt[i] = e
    end

    for (i,e)∈enumerate(eltu)
        new_e = Vector{Int}(undef, 3)
        for j=1:3
            if haskey(interface_dict, vtxu[e[j]])
                new_e[j] = interface_dict[vtxu[e[j]]]
            else
                new_e[j] = e[j] + offset
            end
        end
        elt[i+length(eltl)] = Tuple(new_e)
    end

    return vtx, elt

end
