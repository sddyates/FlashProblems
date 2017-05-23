!!****if* source/Simulation/SimulationMain/ParkerWind/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!!***

!!REORDER(4): solnData, face[xy]Data

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_putPointData,      &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Input Arguments --------------
  integer, intent(in) :: blockID
  !!$ ------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: axis
  real :: rhoZone, presZone
  real :: velxZone, velyZone, velzZone 
  real :: enerZone, ekinZone, eintZone
  real, allocatable, dimension(:) :: xCoord,yCoord,zCoord

  ! Constants.
  real :: pi, kb, mp, Big_G, Msun, Rsun, Mj, Rj
  
  ! Geometric quantities.
  real :: x, y, z, rs, thetas, phis

  ! Stellar parameters.
  real :: Tstar, Mstar, Rstar, RHstar, Pstar

  ! Orbital parameters.
  real :: a, omega_fr, omegas, omega_orb

  ! Stellar wind parameters.
  real :: v_escs, css, rcs, v_ws, v_ws_x, v_ws_y, v_ws_z
  real :: RHOs, PRs, PRSs 

  ! Variables used in the Newton-Raphson routine.
  real :: lambda, eta, b, psi, step, f, df
  integer :: o

  ! Simulation units.
  real :: UNIT_DENSITY, UNIT_LENGTH, UNIT_VELOCITY
  real :: UNIT_MP, UNIT_G, UNIT_TIME, UNIT_MASS

  ! Assigning constants.
  pi = 3.1415926536
  kb = 1.38064852e-16 ![erg K-1]
  mp = 1.67262171e-24 ![g]
  Big_G = 6.67259e-8  ![cm3 g-1 s-2]
  Msun = 1.989e+33    ![g]
  Rsun = 6.955e+10    ![cm]
  Mj = 0.0009543*Msun ![g]
  Rj = 0.10045*Rsun   ![cm]

  ! Defining simulation units.
  UNIT_DENSITY = 1.0e-15
  UNIT_LENGTH = 6.955e+10
  UNIT_VELOCITY = 1.0e+5
  UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY
  UNIT_MASS = UNIT_DENSITY*UNIT_LENGTH**3
  UNIT_G = Big_G*UNIT_DENSITY*UNIT_TIME**2

  ! Normalising simulation parameters.
  a = sim_sep
  Mstar = sim_Ms*Msun/UNIT_MASS
  Rstar = sim_Rs*Rsun/UNIT_LENGTH
  RHstar = sim_RHs/UNIT_DENSITY
  Tstar = sim_Ts
  v_escs = sqrt(2.0*UNIT_G*Mstar/Rstar)
  css = sqrt((2.0*kb*Tstar)/mp)/UNIT_VELOCITY
  Pstar = css*css*RHstar/sim_gamma
  rcs = UNIT_G*Mstar/(2.0*css*css)

  omega_orb = sqrt(UNIT_G*Mstar/(a**3))
  omegas    = omega_orb      
  omega_fr  = omega_orb

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH, IAXIS) - blkLimitsGC(LOW, IAXIS) + 1
  sizeY = blkLimitsGC(HIGH, JAXIS) - blkLimitsGC(LOW, JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH, KAXIS) - blkLimitsGC(LOW, KAXIS) + 1

  allocate(xCoord(sizeX), stat=istat)
  allocate(yCoord(sizeY), stat=istat)
  allocate(zCoord(sizeZ), stat=istat)

  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockID, CENTER, sim_gCell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockID, CENTER, sim_gCell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockID, CENTER, sim_gCell, xCoord, sizeX)
  !------------------------------------------------------------------------------

  ! Loop over cells in the block.
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
    z = zCoord(k)
    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      y = yCoord(j)
      do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        x = xCoord(i)

        rs = sqrt(x**2 + y**2 + z**2)
        thetas = acos(z/rs);
        phis   = atan2(y,x);

        !rp = sqrt((x - a)**2 + y**2 + z**2)
        !thetap = acos(z/rp);
        !phip = atan2(y,(x-a));

        if (rs < 0.5*Rstar) then
          rhoZone = RHstar
          presZone = Pstar + (2.0/3.0)*pi*UNIT_G*RHstar*RHstar*(Rstar*Rstar - rs*rs)
          velxZone = 0.0
          velyZone = 0.0
          velzZone = 0.0
        else if (rs >= 0.5*Rstar .and. rs < Rstar) then
          rhoZone = RHstar
          presZone = Pstar + (2.0/3.0)*pi*UNIT_G*RHstar*RHstar*(Rstar*Rstar - rs*rs)
          velxZone = 0.0
          velyZone = 0.0
          velzZone = 0.0
        else

          ! Newton raphson method for the stellar wind.
          ! This code solves Parkers equation for an 
          ! isothermal stellar wind (Parker 1958).
          lambda = 0.5*(v_escs/css)**2
          eta = rs/Rstar
          b = -3.0 - 4.0*log(lambda/2.0) + 4.0*log(eta)+2.0*(lambda/eta)
          o = 0
          step = 1.0

          if (rs <= rcs) then
            psi = 2.e-8
          else
            psi = 2.5
          endif

          do while (abs(step) > 1.0e-8 .and. o < 1000000) 
            f = -psi + log(psi) + b
            df = -1.0 + (1.0/psi)
            step = f/df
            psi = psi - step
            o = o + 1
          end do

          v_ws = css*sqrt(psi)
          if (isnan(v_ws)) then
            v_ws = 0.0
          endif

          velxZone = sin(thetas)*(v_ws*cos(phis) + sin(phis)*rs*(omega_fr + omegas))
          velyZone = sin(thetas)*(v_ws*sin(phis) - cos(phis)*rs*(omega_fr + omegas))
          velzZone = v_ws*cos(thetas)
          presZone = Pstar*exp(lambda*(Rstar/rs - 1.0) - 0.5*(v_ws/css)**2)
          rhoZone = (RHstar/Pstar)*presZone

        endif

        ekinZone = 0.5*(velxZone*velxZone + velyZone*velyZone + velzZone*velzZone)
        eintZone = presZone/(rhoZone*(sim_gamma - 1.0))
        enerZone = max(eintZone + ekinZone, sim_smallP)

        axis(IAXIS) = i
        axis(JAXIS) = j
        axis(KAXIS) = k

        call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
        call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
        call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
        call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
        call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)
        call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)
        call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
        call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)

      enddo
    enddo
  enddo

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock
