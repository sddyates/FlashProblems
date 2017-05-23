!!****f* source/Simulation/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklist - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)


#include "constants.h"
#include "Flash.h"

  use Simulation_data

  use Grid_interface, ONLY: Grid_getCellCoords
  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_getBlkPtr
  use Grid_interface, ONLY: Grid_releaseBlkPtr
  use Grid_interface, ONLY: Grid_getBlkBC

  implicit none

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real,    intent(in) :: dt
  real,    intent(in) :: stime
  integer :: lb

  real, pointer, dimension(:,:,:,:) :: blkPtr!, facexData, faceyData, facezData

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: rhoZone, presZone
  real :: velxZone, velyZone, velzZone 
  real :: magxZone, magyZone, magzZone
  real :: magpZone, divbZone
  real :: enerZone, ekinZone, eintZone
  real, allocatable, dimension(:) :: xCoord,yCoord,zCoord

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
  omegas    = omega_orb;      
  omega_fr  = omega_orb;

  do lb = 1, blkcnt
    call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blklst(lb), blkPtr)
 
    sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
    sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
    sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
 
    allocate(xCoord(sizeX),stat=istat)
    allocate(yCoord(sizeY),stat=istat)
    allocate(zCoord(sizeZ),stat=istat)
 
    xCoord = 0.0
    yCoord = 0.0
    zCoord = 0.0

    if (NDIM == 3) call Grid_getCellCoords(KAXIS,blklst(lb),CENTER,sim_gCell,zCoord,sizeZ)
    if (NDIM >= 2) call Grid_getCellCoords(JAXIS,blklst(lb),CENTER,sim_gCell,yCoord,sizeY)
    call Grid_getCellCoords(IAXIS,blklst(lb),CENTER,sim_gCell,xCoord,sizeX)


!#if NFACE_VARS > 0
  !if (sim_killdivb) then
!     call Grid_getBlkPtr(blklst(lb),facexData,FACEX)
!     call Grid_getBlkPtr(blklst(lb),faceyData,FACEY)
!     call Grid_getBlkPtr(blklst(lb),facezData,FACEZ)
  !endif
!#endif


    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      z = zCoord(k)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        y = yCoord(j)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          x = xCoord(i)

          rs = sqrt(x**2 + y**2 + z**2)
          thetas = acos(z/rs);
          phis   = atan2(y,x);

          !rp = sqrt((x - a)**2 + y**2 + z**2)
          !thetap = acos(z/rp);
          !phip = atan2(y,(x-a));

          if (rs < 0.5*Rstar) then
            blkPtr(DENS_VAR,i,j,k) = RHstar
            blkPtr(PRES_VAR,i,j,k) = Pstar + (2.0/3.0)*pi*UNIT_G*RHstar*RHstar*(Rstar*Rstar - rs*rs)
            blkPtr(VELX_VAR,i,j,k) = 0.0
            blkPtr(VELY_VAR,i,j,k) = 0.0
            blkPtr(VELZ_VAR,i,j,k) = 0.0
            blkPtr(MAGX_VAR,i,j,k) = 0.0
            blkPtr(MAGY_VAR,i,j,k) = 0.0
            blkPtr(MAGZ_VAR,i,j,k) = 16.0*B0star
!#if NFACE_VARS > 0
!            if (sim_killdivb) then
!              facexData(MAG_FACE_VAR,i,j,k) = 0.0
!              faceyData(MAG_FACE_VAR,i,j,k) = 0.0
!              facezData(MAG_FACE_VAR,i,j,k) = 16.0*B0star
!            endif
!#endif
            magpZone = 0.5*(blkPtr(MAGX_VAR,i,j,k)*blkPtr(MAGX_VAR,i,j,k) & 
                          + blkPtr(MAGY_VAR,i,j,k)*blkPtr(MAGY_VAR,i,j,k) &
                          + blkPtr(MAGZ_VAR,i,j,k)*blkPtr(MAGZ_VAR,i,j,k))
            ekinZone = 0.5*(blkPtr(VELX_VAR,i,j,k)*blkPtr(VELX_VAR,i,j,k) &
                          + blkPtr(VELY_VAR,i,j,k)*blkPtr(VELY_VAR,i,j,k) & 
                          + blkPtr(VELZ_VAR,i,j,k)*blkPtr(VELZ_VAR,i,j,k))
            eintZone = blkPtr(PRES_VAR,i,j,k)/(blkPtr(DENS_VAR,i,j,k)*(sim_gamma - 1.0))
            enerZone = max(eintZone + ekinZone, sim_smallP)
            divbZone = 0.0

            blkPtr(DIVB_VAR,i,j,k) = divbZone 
            blkPtr(MAGP_VAR,i,j,k) = magpZone
            blkPtr(ENER_VAR,i,j,k) = enerZone 
            blkPtr(EINT_VAR,i,j,k) = eintZone
            blkPtr(GAMC_VAR,i,j,k) = sim_gamma
            blkPtr(GAME_VAR,i,j,k) = sim_gamma
          else if (rs >= 0.5*Rstar .and. rs < Rstar) then
            blkPtr(DENS_VAR,i,j,k) = RHstar
            blkPtr(PRES_VAR,i,j,k) = Pstar + (2.0/3.0)*pi*UNIT_G*RHstar*RHstar*(Rstar*Rstar - rs*rs)
            blkPtr(VELX_VAR,i,j,k) = 0.0
            blkPtr(VELY_VAR,i,j,k) = 0.0
            blkPtr(VELZ_VAR,i,j,k) = 0.0
            blkPtr(MAGX_VAR,i,j,k) = 3.0*x*z*B0star*Rstar**3*rs**(-5)
            blkPtr(MAGY_VAR,i,j,k) = 3.0*z*y*B0star*Rstar**3*rs**(-5)
            blkPtr(MAGZ_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0star*Rstar**3*rs**(-5)
!#if NFACE_VARS > 0
!            if (sim_killdivb) then
!              facexData(MAG_FACE_VAR,i,j,k) = 3.0*x*z*B0star*Rstar**3*rs**(-5) 
!              faceyData(MAG_FACE_VAR,i,j,k) = 3.0*z*y*B0star*Rstar**3*rs**(-5)
!              facezData(MAG_FACE_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0star*Rstar**3*rs**(-5)
!            endif
!#endif
            magpZone = 0.5*(blkPtr(MAGX_VAR,i,j,k)*blkPtr(MAGX_VAR,i,j,k) & 
                          + blkPtr(MAGY_VAR,i,j,k)*blkPtr(MAGY_VAR,i,j,k) &
                          + blkPtr(MAGZ_VAR,i,j,k)*blkPtr(MAGZ_VAR,i,j,k))
            ekinZone = 0.5*(blkPtr(VELX_VAR,i,j,k)*blkPtr(VELX_VAR,i,j,k) &
                          + blkPtr(VELY_VAR,i,j,k)*blkPtr(VELY_VAR,i,j,k) & 
                          + blkPtr(VELZ_VAR,i,j,k)*blkPtr(VELZ_VAR,i,j,k))
            eintZone = blkPtr(PRES_VAR,i,j,k)/(blkPtr(DENS_VAR,i,j,k)*(sim_gamma - 1.0))
            enerZone = max(eintZone + ekinZone, sim_smallP)
            divbZone = 0.0
            blkPtr(DIVB_VAR,i,j,k) = divbZone 
            blkPtr(MAGP_VAR,i,j,k) = magpZone
            blkPtr(ENER_VAR,i,j,k) = enerZone 
            blkPtr(EINT_VAR,i,j,k) = eintZone
            blkPtr(GAMC_VAR,i,j,k) = sim_gamma
            blkPtr(GAME_VAR,i,j,k) = sim_gamma
          else if (rs >= Rstar .and. rs < 1.5*Rstar) then

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

            blkPtr(VELX_VAR,i,j,k) = sin(thetas)*(v_ws*cos(phis) + sin(phis)*rs*(omega_fr + omegas))
            blkPtr(VELY_VAR,i,j,k) = sin(thetas)*(v_ws*sin(phis) - cos(phis)*rs*(omega_fr + omegas))
            blkPtr(VELZ_VAR,i,j,k) = v_ws*cos(thetas)
            blkPtr(PRES_VAR,i,j,k) = Pstar*exp(lambda*(Rstar/rs - 1.0) - 0.5*(v_ws/css)**2)
            blkPtr(DENS_VAR,i,j,k) = (RHstar/Pstar)*blkPtr(PRES_VAR,i,j,k)

            ekinZone = 0.5*(blkPtr(VELX_VAR,i,j,k)*blkPtr(VELX_VAR,i,j,k) &
                          + blkPtr(VELY_VAR,i,j,k)*blkPtr(VELY_VAR,i,j,k) & 
                          + blkPtr(VELZ_VAR,i,j,k)*blkPtr(VELZ_VAR,i,j,k))
            eintZone = blkPtr(PRES_VAR,i,j,k)/(blkPtr(DENS_VAR,i,j,k)*(sim_gamma - 1.0))
            enerZone = max(eintZone + ekinZone, sim_smallP)
            divbZone = 0.0

            blkPtr(DIVB_VAR,i,j,k) = divbZone 
            blkPtr(ENER_VAR,i,j,k) = enerZone 
            blkPtr(EINT_VAR,i,j,k) = eintZone
            blkPtr(GAMC_VAR,i,j,k) = sim_gamma
            blkPtr(GAME_VAR,i,j,k) = sim_gamma
          endif

        end do
      end do
    end do

  end do

  call Grid_releaseBlkPtr(blklst(lb), blkPtr)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_adjustEvolution
