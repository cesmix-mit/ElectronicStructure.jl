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
    linearize(energies::Vector{Unitful.Energy})
    
Linearize a vector of energies. Units are removed.
"""
function linearize(energies::Vector{Unitful.Energy})
    return map(x->x.val, energies)
end

"""
    linearize(forces::Vector{Vector{SVector{D, Unitful.Force}}}) where {D}
    
Linearize a vector of forces. Units are removed.
"""
function linearize_aux(v::Vector)
    return Vector(map(x->x.val, vcat(v...)))
end
function linearize(forces::Vector{Vector{SVector{D, Unitful.Force}}}) where {D}
    return vcat(linearize_aux.(forces)...)
end

"""
    linearize(stresses::Vector{SMatrix{D, D, Unitful.Pressure}}) where {D}
    
Linearize a vector of stresses. Only the upper triangular part of each stress is
considered. Units are removed.
"""
function linearize(stresses::Vector{SMatrix{D, D, Unitful.Pressure}}) where {D}
    return map(x->x.val, vcat(map(vec, lower_tri_vec.(stresses))...))
end

"""
    linearize(d::SmallESData)
    
Linearize the energies, forces, and stresses of a SmallESData struct. 
Only the upper triangular part of each stress is considered. Units are removed.
"""
function linearize(d::SmallESData)
    return vcat( linearize(d.energies),
                 linearize(d.forces),
                 linearize(d.stresses))
end

