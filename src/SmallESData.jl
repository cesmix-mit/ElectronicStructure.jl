export SmallESData, gen_test_data, linearize

struct SmallESData{D} <: ElectronicStructureData{D}
    energies::Vector{Unitful.Energy}
    forces::Vector{Vector{SVector{D, Unitful.Force}}}
    stresses::Vector{SMatrix{D, D, Unitful.Pressure}} # TODO: stress type in Unitful?
end

"""
    gen_test_data(D, atomic_confs, p)

Generate test DFT data
"""
function gen_test_data(D, atomic_confs, p)
    pe  = [ potential_energy(s, p) for s in atomic_confs ]
    f   = [ forces(s, p) for s in atomic_confs ]
    s   = [ virial_stress(s, p) for s in atomic_confs ]
    data = SmallESData{D}(pe, f, s)
    return data
end


"""
    lower_tri_vec(s::SMatrix{D, D}) where {D}
    
 Returns a vector containing the lower triangle of a matrix.
"""
function lower_tri_vec(s::SMatrix{D, D}) where {D}
    return [ s[i, j]  for j in 1:3 for i in j:3 ]
end

"""
    linearize(v::Vector)
    
Linearize a vector of forces, or a vector of stresses. Units are removed.
"""
function linearize(v::Vector)
    return map(x->x.val, vcat(map(vec, v)...))
end

"""
    linearize(d::SmallESData)
    
Linearize the energies, forces, and stresses of a SmallESData struct. Units are removed.
"""
function linearize(d::SmallESData)
    return vcat( map(x->x.val, d.energies), 
                 vcat(linearize.(d.forces)...),
                 linearize(lower_tri_vec.(d.stresses)))
end



