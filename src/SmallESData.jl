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
    linearize(v::Vector)
    
Linearize a vector of forces, or a vector of stresses. Units are removed.
"""
function linearize(v::Vector)
    return map(x->x.val, reduce(vcat, map(vec, v)))
end

