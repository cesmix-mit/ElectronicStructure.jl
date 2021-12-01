export SmallESData, gen_test_data

struct SmallESData{D} <: ElectronicStructureData{D}
    energies::Vector{Unitful.Energy}
    forces::Vector{Vector{SVector{D, Unitful.Force}}}
    stresses::Vector{SMatrix{D, D, Float64}} # Q: stress type in Unitful?
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

