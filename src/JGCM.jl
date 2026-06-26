module JGCM

# ============================================================================== #
# CORE STRUCTURES
# ============================================================================== #
include("Core/Gauss_And_Legendre.jl")
include("Core/Spectral_Spherical_Mesh.jl")
include("Core/Atmo_Data.jl")
include("Core/Vert_Coordinate.jl")
include("Core/Dyn_Data.jl")
include("Core/Variable_Mappings.jl")
include("Core/Experiment_Configuration.jl")

# ============================================================================== #
# DYNAMICS & PHYSICS
# ============================================================================== #
# Time Integration
include("Dynamics/Time_Integrator.jl")

# Dynamics Components
include("Dynamics/Press_And_Geopot.jl")
include("Dynamics/Semi_Implicit.jl")
include("Dynamics/Barotropic_Dynamics.jl")
include("Dynamics/Shallow_Water_Dynamics.jl")

# Physics Parameterizations
include("Physics/Atmos_Param.jl")

# Full Spectral Dynamics (Depends on all above)
include("Dynamics/Spectral_Dynamics.jl")

# ============================================================================== #
# INITIALIZATION & OUTPUT
# ============================================================================== #
include("Output/Restart_Manager.jl")
include("Initialization/Initial_Conditions.jl")
include("Output/Output_Mappings.jl")
include("Output/Vertical_Interpolation.jl")
include("Output/Output_Manager.jl")

# ============================================================================== #
# DRIVER
# ============================================================================== #
include("Driver.jl")

# ============================================================================== #
# EXPORTS
# ============================================================================== #
# Core Types
using .Gauss_And_Legendre_Module
export Compute_Gaussian, Compute_Legendre

using .Spectral_Spherical_Mesh_Module
export Spectral_Spherical_Mesh
export Trans_Spherical_To_Grid!, Trans_Grid_To_Spherical!, Trans_Grid_To_Fourier!
export Vor_Div_From_Grid_UV!, UV_Grid_From_Vor_Div!
export Compute_Gradient_Cos!, Add_Horizontal_Advection!
export Divide_By_Cos!, Multiply_By_Cos!

using .Atmo_Data_Module
export Atmo_Data

using .Vert_Coordinate_Module
export Vert_Coordinate

using .Dyn_Data_Module
export Dyn_Data

using .Variable_Mappings_Module
export Get_Data_Pointer

using .Experiment_Configuration
export Model_Config

# Initialization
using .Initial_Conditions
export Initialize_Atmos_State!

# Dynamics (Exporting explicitly for testing/access)
using .Spectral_Dynamics_Module
export Spectral_Initialize_Fields!, Get_Topography!, Spectral_Dynamics!

using .Barotropic_Dynamics_Module
export Barotropic_Dynamics!

using .Shallow_Water_Dynamics_Module
export Shallow_Water_Dynamics!

# Output
using .Restart_Manager_Module
export Restart_Manager, Write_Restart_File, Load_Restart_File!
using .Output_Manager_Module
export Output_Manager, Update_Output!, Finalize_Output!

# Main Driver
using .Driver
export JGCM_Simulate

end
