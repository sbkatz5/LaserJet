!!****if* source/Simulation/SimulationMain/LaserSlab/Simulation_init
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
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for a particular simulation
!!
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init()
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  
  implicit none

#include "constants.h"
#include "Flash.h"

  real :: xmin, xmax, ymin, ymax
  integer :: lrefine_max, nblockx, nblocky
  character(len=MAX_STRING_LENGTH) :: str

  call RuntimeParameters_get('sim_windowRadius', sim_windowRadius)
  call RuntimeParameters_get('sim_windowHeight', sim_windowHeight)

  call RuntimeParameters_get('sim_windowWidth', sim_windowWidth)
  call RuntimeParameters_get('sim_windowRad', sim_windowRad)
  call RuntimeParameters_get('sim_rhoWind', sim_rhoWind)
  call RuntimeParameters_get('sim_teleWind', sim_teleWind)
  call RuntimeParameters_get('sim_tionWind', sim_tionWind)
  call RuntimeParameters_get('sim_tradWind', sim_tradWind)

  call RuntimeParameters_get('sim_LaserEntry', sim_LaserEntry)
  call RuntimeParameters_get('sim_LaserExit', sim_LaserExit)
  call RuntimeParameters_get('sim_rhoVacu', sim_rhoVacu)
  call RuntimeParameters_get('sim_teleVacu', sim_teleVacu)
  call RuntimeParameters_get('sim_tionVacu', sim_tionVacu)
  call RuntimeParameters_get('sim_tradVacu', sim_tradVacu)

  call RuntimeParameters_get('sim_washerWidth', sim_washerWidth)
  call RuntimeParameters_get('sim_washerOuterRad', sim_washerOuterRad)
  call RuntimeParameters_get('sim_washerInnerRad', sim_washerInnerRad)
  call RuntimeParameters_get('sim_rhoWash', sim_rhoWash)
  call RuntimeParameters_get('sim_teleWash', sim_teleWash)
  call RuntimeParameters_get('sim_tionWash', sim_tionWash)
  call RuntimeParameters_get('sim_tradWash', sim_tradWash)

  call RuntimeParameters_get('smallX', sim_smallX)

  call RuntimeParameters_get('sim_initGeom', sim_initGeom)

#ifdef FLASH_USM_MHD
  call RuntimeParameters_get('killdivb', sim_killdivb)
#endif
end subroutine Simulation_init
