################################################################################
#
#    Interface.jl
#
################################################################################

using Unitful, StaticArrays

export ElectronicStructureData, gen_test_data

# Abstract types ###############################################################
abstract type ElectronicStructureData{D} end

# Functions ####################################################################
function gen_test_data end
