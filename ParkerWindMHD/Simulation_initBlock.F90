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
  real :: magxZone, magyZone, magzZone
  real :: magpZone, divbZone
  real :: enerZone, ekinZone, eintZone
  real, allocatable, dimension(:) :: xCoord,yCoord,zCoord
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData, facezData

  ! Constants.
  real :: pi, kb, mp, Big_G, Msun, Rsun, Mj, Rj
  
  ! Geometric quantities.
  real :: x, y, z, rs, thetas, phis

  ! Stellar parameters.
  real :: Tstar, Mstar, Rstar, RHstar, Pstar, B0star

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
  real :: UNIT_MP, UNIT_G, UNIT_TIME, UNIT_MASS, UNIT_MAG

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
  UNIT_MAG = sqrt(4.0*pi*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)
  UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY
  UNIT_MASS = UNIT_DENSITY*UNIT_LENGTH**3
  UNIT_G = Big_G*UNIT_DENSITY*UNIT_TIME**2

  ! Normalising simulation parameters.
  a = sim_sep
  Mstar = sim_Ms*Msun/UNIT_MASS
  Rstar = sim_Rs*Rsun/UNIT_LENGTH
  RHstar = sim_RHs/UNIT_DENSITY
  Tstar = sim_Ts
  B0star = sim_B0s/UNIT_MAG
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

  allocate(xCoord(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))

  xCoord  = 0.0
  yCoord  = 0.0
  zCoord  = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockID, CENTER, sim_gCell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockID, CENTER, sim_gCell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockID, CENTER, sim_gCell, xCoord, sizeX)
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  !if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
  !endif
#endif

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
          solnData(DENS_VAR,i,j,k) = RHstar
          solnData(PRES_VAR,i,j,k) = Pstar + (2.0/3.0)*pi*UNIT_G*RHstar*RHstar*(Rstar*Rstar - rs*rs)
          solnData(VELX_VAR,i,j,k) = 0.0
          solnData(VELY_VAR,i,j,k) = 0.0
          solnData(VELZ_VAR,i,j,k) = 0.0
          solnData(MAGX_VAR,i,j,k) = 0.0
          solnData(MAGY_VAR,i,j,k) = 0.0
          solnData(MAGZ_VAR,i,j,k) = 16.0*B0star
#if NFACE_VARS > 0
          if (sim_killdivb) then
            facexData(MAG_FACE_VAR,i,j,k) = 0.0
            faceyData(MAG_FACE_VAR,i,j,k) = 0.0
            facezData(MAG_FACE_VAR,i,j,k) = 16.0*B0star
          endif
#endif
        else if (rs >= 0.5*Rstar .and. rs < Rstar) then
          solnData(DENS_VAR,i,j,k) = RHstar
          solnData(PRES_VAR,i,j,k) = Pstar + (2.0/3.0)*pi*UNIT_G*RHstar*RHstar*(Rstar*Rstar - rs*rs)
          solnData(VELX_VAR,i,j,k) = 0.0
          solnData(VELY_VAR,i,j,k) = 0.0
          solnData(VELZ_VAR,i,j,k) = 0.0
          solnData(MAGX_VAR,i,j,k) = 3.0*x*z*B0star*Rstar**3*rs**(-5) 
          solnData(MAGY_VAR,i,j,k) = 3.0*z*y*B0star*Rstar**3*rs**(-5) 
          solnData(MAGZ_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0star*Rstar**3*rs**(-5) 

#if NFACE_VARS > 0
          if (sim_killdivb) then
            facexData(MAG_FACE_VAR,i,j,k) = 3.0*x*z*B0star*Rstar**3*rs**(-5) 
            faceyData(MAG_FACE_VAR,i,j,k) = 3.0*z*y*B0star*Rstar**3*rs**(-5)
            facezData(MAG_FACE_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0star*Rstar**3*rs**(-5)
          endif
#endif
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

          solnData(VELX_VAR,i,j,k) = sin(thetas)*(v_ws*cos(phis) + sin(phis)*rs*(omega_fr + omegas))
          solnData(VELY_VAR,i,j,k) = sin(thetas)*(v_ws*sin(phis) - cos(phis)*rs*(omega_fr + omegas))
          solnData(VELZ_VAR,i,j,k) = v_ws*cos(thetas)
          solnData(PRES_VAR,i,j,k) = Pstar*exp(lambda*(Rstar/rs - 1.0) - 0.5*(v_ws/css)**2)
          solnData(DENS_VAR,i,j,k) = (RHstar/Pstar)*solnData(PRES_VAR,i,j,k)

          solnData(MAGX_VAR,i,j,k) = 3.0*x*z*B0star*Rstar**3*rs**(-5)
          solnData(MAGY_VAR,i,j,k) = 3.0*z*y*B0star*Rstar**3*rs**(-5)
          solnData(MAGZ_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0star*Rstar**3*rs**(-5)

#if NFACE_VARS > 0
          if (sim_killdivb) then
            facexData(MAG_FACE_VAR,i,j,k) = 3.0*x*z*B0star*Rstar**3*rs**(-5) 
            faceyData(MAG_FACE_VAR,i,j,k) = 3.0*z*y*B0star*Rstar**3*rs**(-5)
            facezData(MAG_FACE_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0star*Rstar**3*rs**(-5)
          endif
#endif

        endif

        magpZone = 0.5*(solnData(MAGX_VAR,i,j,k)*solnData(MAGX_VAR,i,j,k) &
                      + solnData(VELY_VAR,i,j,k)*solnData(VELY_VAR,i,j,k) &
                      + solnData(MAGZ_VAR,i,j,k)*solnData(MAGZ_VAR,i,j,k))
        ekinZone = 0.5*(solnData(VELX_VAR,i,j,k)*solnData(VELX_VAR,i,j,k) &
                      + solnData(VELY_VAR,i,j,k)*solnData(VELY_VAR,i,j,k) &
                      + solnData(VELZ_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k))
        eintZone = solnData(PRES_VAR,i,j,k)/(solnData(DENS_VAR,i,j,k)*(sim_gamma - 1.0))
        enerZone = max(eintZone + ekinZone, sim_smallP)
        divbZone = 0.0

        solnData(DIVB_VAR,i,j,k) = divbZone
        solnData(MAGP_VAR,i,j,k) = magpZone
        solnData(ENER_VAR,i,j,k) = enerZone
        solnData(GAME_VAR,i,j,k) = sim_gamma
        solnData(GAMC_VAR,i,j,k) = sim_gamma

        axis(IAXIS) = i
        axis(JAXIS) = j
        axis(KAXIS) = k

        !call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
        !call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)

        !call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
        !call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
        !call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

        !call Grid_putPointData(blockId, CENTER, MAGX_VAR, EXTERIOR, axis, magxZone)
        !call Grid_putPointData(blockId, CENTER, MAGY_VAR, EXTERIOR, axis, magyZone)
        !call Grid_putPointData(blockId, CENTER, MAGZ_VAR, EXTERIOR, axis, magzZone)

        !call Grid_putPointData(blockId, CENTER, DIVB_VAR, EXTERIOR, axis, divbZone)
        !call Grid_putPointData(blockId, CENTER, MAGP_VAR, EXTERIOR, axis, magpZone)
        !call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)
        !call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
        !call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)

      enddo
    enddo
  enddo

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock
