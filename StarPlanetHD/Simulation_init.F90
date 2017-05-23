!!****if* source/Simulation/SimulationMain/ParkerWind/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! ARGUMENTS
!!
!!   
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!!  Initializes initial conditions for Parker Wind simulation.
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  call RuntimeParameters_get('gamma', sim_gamma)

  call RuntimeParameters_get('seperation', sim_sep)
  call RuntimeParameters_get('star_temp', sim_Ts)
  call RuntimeParameters_get('star_mass', sim_Ms)
  call RuntimeParameters_get('star_radius', sim_Rs)
  call RuntimeParameters_get('star_mag', sim_B0s)
  call RuntimeParameters_get('star_rho', sim_RHs)

  call RuntimeParameters_get('planet_temp', sim_Tp)
  call RuntimeParameters_get('planet_mass', sim_Mp)
  call RuntimeParameters_get('planet_radius', sim_Rp)
  call RuntimeParameters_get('planet_mag', sim_B0p)
  call RuntimeParameters_get('planet_rho', sim_RHp)
    
  call RuntimeParameters_get('xmin', sim_xMin)
  call RuntimeParameters_get('ymin', sim_yMin)
  call RuntimeParameters_get('zmin', sim_zMin)

  call RuntimeParameters_get('xmax', sim_xMax)
  call RuntimeParameters_get('ymax', sim_yMax)
  call RuntimeParameters_get('zmax', sim_zMax)

  call RuntimeParameters_get('xCtr', sim_xCtr)
  call RuntimeParameters_get('yCtr', sim_yCtr)
  call RuntimeParameters_get('zCtr', sim_zCtr)

  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX)

  sim_gCell = .true.
  
end subroutine Simulation_init
