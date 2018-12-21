!!****if* source/Simulation/SimulationMain/LaserSlab/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!! 
!!***
module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!  
  real, save :: sim_windowRadius
  real, save :: sim_windowHeight

  real,    save :: sim_windowWidth
  real,    save :: sim_windowRad
  real,    save :: sim_rhoWind  
  real,    save :: sim_teleWind 
  real,    save :: sim_tionWind 
  real,    save :: sim_tradWind 
  real,    save :: sim_zminWind
  integer, save :: sim_eosWind

  real,    save :: sim_LaserEntry
  real,    save :: sim_LaserExit
  real,    save :: sim_rhoVacu  
  real,    save :: sim_teleVacu 
  real,    save :: sim_tionVacu 
  real,    save :: sim_tradVacu 
  integer, save :: sim_eosVacu

  real,    save :: sim_washerWidth
  real,    save :: sim_washerOuterRad
  real,    save :: sim_washerInnerRad
  real,    save :: sim_rhoWash
  real,    save :: sim_teleWash
  real,    save :: sim_tionWash
  real,    save :: sim_tradWash
  real,    save :: sim_zminWash
  integer, save :: sim_eosWash

  logical, save :: sim_killdivb = .FALSE.
  real, save :: sim_smallX
  character(len=MAX_STRING_LENGTH), save :: sim_initGeom


end module Simulation_data


