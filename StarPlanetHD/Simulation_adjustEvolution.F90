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

  real, pointer :: blkPtr(:,:,:,:)

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: rhoZone, presZone
  real :: velxZone, velyZone, velzZone 
  real :: enerZone, ekinZone, eintZone
  real, allocatable, dimension(:) :: xCoord,yCoord,zCoord

  ! Constants.
  real :: pi, kb, mp, Big_G, Msun, Rsun, Mj, Rj
  
  ! Geometric quantities.
  real :: x, y, z, rs, thetas, phis, rp, thetap, phip

  ! Stellar parameters.
  real :: Tstar, Mstar, Rstar, RHstar, Pstar

  ! Planetary parameters.
  real :: Tplanet, Mplanet, Rplanet, RHplanet, Pplanet

  ! Orbital parameters.
  real :: a, omega_fr, omegas, omegap, omega_orb

  ! Stellar wind parameters.
  real :: v_escs, css, rcs, v_ws, v_ws_x, v_ws_y, v_ws_z
  real :: RHOs, PRs, PRSs 

  ! Planetary wind parameters.
  real :: v_escp, csp, rcp, v_wp, v_wp_x, v_wp_y, v_wp_z
  real :: RHOp, PRp, PRSp

  ! Variables used in the Newton-Raphson routine.
  real :: lambda, eta, b, psi, step, f, df
  real :: lambdas, etas
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
  a = sim_sep*0.01*1.496e+13/UNIT_LENGTH
  Mstar = sim_Ms*Msun/UNIT_MASS
  Rstar = sim_Rs*Rsun/UNIT_LENGTH
  RHstar = sim_RHs/UNIT_DENSITY
  Tstar = sim_Ts
  v_escs = sqrt(2.0*UNIT_G*Mstar/Rstar)
  css = sqrt((2.0*kb*Tstar)/mp)/UNIT_VELOCITY
  Pstar = css*css*RHstar/sim_gamma
  rcs = UNIT_G*Mstar/(2.0*css*css)

  Mplanet = sim_Mp*Mj/UNIT_MASS
  Rplanet = sim_Rp*Rj/UNIT_LENGTH
  RHplanet = sim_RHp/UNIT_DENSITY
  Tplanet = sim_Tp
  v_escp = sqrt(2.0*UNIT_G*Mplanet/Rplanet)
  csp = sqrt((2.0*kb*Tplanet)/mp)/UNIT_VELOCITY
  Pplanet = csp*csp*RHplanet/sim_gamma
  rcp = UNIT_G*Mplanet/(2.0*csp*csp)

  omega_orb = sqrt(UNIT_G*Mstar/(a**3))
  omegas    = omega_orb      
  omegap    = omega_orb      
  omega_fr  = omega_orb

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

    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      z = zCoord(k)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        y = yCoord(j)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          x = xCoord(i)

          rs = sqrt(x**2 + y**2 + z**2)
          thetas = acos(z/rs);
          phis   = atan2(y,x);

          rp = sqrt((x - a)**2 + y**2 + z**2)
          thetap = acos(z/rp);
          phip = atan2(y,(x-a));

          if (rs < 0.5*Rstar) then
            rhoZone = RHstar
            presZone = Pstar + (2.0/3.0)*pi*UNIT_G*RHstar*RHstar*(Rstar*Rstar - rs*rs)
            velxZone = 0.0
            velyZone = 0.0
            velzZone = 0.0
            ekinZone = 0.5*(velxZone*velxZone + velyZone*velyZone + velzZone*velzZone)
            eintZone = presZone/(rhoZone*(sim_gamma - 1.0))
            enerZone = max(eintZone + ekinZone, sim_smallP)
            blkPtr(DENS_VAR,i,j,k) = rhoZone
            blkPtr(PRES_VAR,i,j,k) = presZone
            blkPtr(VELX_VAR,i,j,k) = velxZone
            blkPtr(VELY_VAR,i,j,k) = velyZone
            blkPtr(VELZ_VAR,i,j,k) = velzZone
            blkPtr(ENER_VAR,i,j,k) = enerZone 
            blkPtr(EINT_VAR,i,j,k) = eintZone
            blkPtr(GAMC_VAR,i,j,k) = sim_gamma
            blkPtr(GAME_VAR,i,j,k) = sim_gamma
          else if (rs >= 0.5*Rstar .and. rs < Rstar) then
            rhoZone = RHstar
            presZone = Pstar + (2.0/3.0)*pi*UNIT_G*RHstar*RHstar*(Rstar*Rstar - rs*rs)
            velxZone = 0.0
            velyZone = 0.0
            velzZone = 0.0
            ekinZone = 0.5*(velxZone*velxZone + velyZone*velyZone + velzZone*velzZone)
            eintZone = presZone/(rhoZone*(sim_gamma - 1.0))
            enerZone = max(eintZone + ekinZone, sim_smallP)
            blkPtr(DENS_VAR,i,j,k) = rhoZone
            blkPtr(PRES_VAR,i,j,k) = presZone
            blkPtr(VELX_VAR,i,j,k) = velxZone
            blkPtr(VELY_VAR,i,j,k) = velyZone
            blkPtr(VELZ_VAR,i,j,k) = velzZone
            blkPtr(ENER_VAR,i,j,k) = enerZone 
            blkPtr(EINT_VAR,i,j,k) = eintZone
            blkPtr(GAMC_VAR,i,j,k) = sim_gamma
            blkPtr(GAME_VAR,i,j,k) = sim_gamma
          else if (rs >= Rstar .and. rs < 1.5*Rstar) then

            ! Newton raphson method for the stellar wind.
            ! This code solves Parkers equation for an 
            ! isothermal stellar wind (Parker 1958).
            lambdas = 0.5*(v_escs/css)**2
            etas = rs/Rstar
            b = -3.0 - 4.0*log(lambdas/2.0) + 4.0*log(etas)+2.0*(lambdas/etas)
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
            presZone = Pstar*exp(lambdas*(Rstar/rs - 1.0) - 0.5*(v_ws/css)**2)
            rhoZone = (RHstar/Pstar)*presZone
            ekinZone = 0.5*(velxZone*velxZone + velyZone*velyZone + velzZone*velzZone)
            eintZone = presZone/(rhoZone*(sim_gamma - 1.0))
            enerZone = max(eintZone + ekinZone, sim_smallP)
            blkPtr(DENS_VAR,i,j,k) = rhoZone
            blkPtr(PRES_VAR,i,j,k) = presZone
            blkPtr(VELX_VAR,i,j,k) = velxZone
            blkPtr(VELY_VAR,i,j,k) = velyZone
            blkPtr(VELZ_VAR,i,j,k) = velzZone
            blkPtr(ENER_VAR,i,j,k) = enerZone 
            blkPtr(EINT_VAR,i,j,k) = eintZone
            blkPtr(GAMC_VAR,i,j,k) = sim_gamma
            blkPtr(GAME_VAR,i,j,k) = sim_gamma
          endif

          if (rp < 0.5*Rplanet) then
            rhoZone = RHplanet
            presZone = Pplanet + (2.0/3.0)*pi*UNIT_G*RHplanet*RHplanet*(Rplanet*Rplanet - rp*rp)
            velxZone = 0.0
            velyZone = 0.0
            velzZone = 0.0
            ekinZone = 0.5*(velxZone*velxZone + velyZone*velyZone + velzZone*velzZone)
            eintZone = presZone/(rhoZone*(sim_gamma - 1.0))
            enerZone = max(eintZone + ekinZone, sim_smallP)
            blkPtr(DENS_VAR,i,j,k) = rhoZone
            blkPtr(PRES_VAR,i,j,k) = presZone
            blkPtr(VELX_VAR,i,j,k) = velxZone
            blkPtr(VELY_VAR,i,j,k) = velyZone
            blkPtr(VELZ_VAR,i,j,k) = velzZone
            blkPtr(ENER_VAR,i,j,k) = enerZone 
            blkPtr(EINT_VAR,i,j,k) = eintZone
            blkPtr(GAMC_VAR,i,j,k) = sim_gamma
            blkPtr(GAME_VAR,i,j,k) = sim_gamma
          else if (rp >= 0.5*Rplanet .and. rp < Rplanet) then
            rhoZone = RHplanet
            presZone = Pplanet + (2.0/3.0)*pi*UNIT_G*RHplanet*RHplanet*(Rplanet*Rplanet - rp*rp)
            velxZone = 0.0
            velyZone = 0.0
            velzZone = 0.0
            ekinZone = 0.5*(velxZone*velxZone + velyZone*velyZone + velzZone*velzZone)
            eintZone = presZone/(rhoZone*(sim_gamma - 1.0))
            enerZone = max(eintZone + ekinZone, sim_smallP)
            blkPtr(DENS_VAR,i,j,k) = rhoZone
            blkPtr(PRES_VAR,i,j,k) = presZone
            blkPtr(VELX_VAR,i,j,k) = velxZone
            blkPtr(VELY_VAR,i,j,k) = velyZone
            blkPtr(VELZ_VAR,i,j,k) = velzZone
            blkPtr(ENER_VAR,i,j,k) = enerZone 
            blkPtr(EINT_VAR,i,j,k) = eintZone
            blkPtr(GAMC_VAR,i,j,k) = sim_gamma
            blkPtr(GAME_VAR,i,j,k) = sim_gamma
          else if (rp >= Rplanet .and. rp <= 1.5*Rplanet) then
            ! Newton raphson method for the stellar wind.
            ! This code solves Parkers equation for an 
            ! isothermal stellar wind (Parker 1958).
            lambda = 0.5*(v_escp/csp)**2
            eta = rp/Rplanet
            b = -3.0 - 4.0*log(lambda/2.0) + 4.0*log(eta)+2.0*(lambda/eta)
            o = 0
            step = 1.0
            if (rp <= rcp) then
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
            v_wp = csp*sqrt(psi)
            if (isnan(v_wp)) then
              v_wp = 0.0
            endif
            velxZone = sin(thetap)*(v_wp*cos(phip) + sin(phis)*rs*omega_fr - sin(phip)*rp*omegap)
            velyZone = sin(thetap)*(v_wp*sin(phip) - cos(phis)*rs*omega_fr + a*omega_orb + cos(phip)*rp*omegap)
            velzZone = v_wp*cos(thetap)
            presZone = Pplanet*exp(lambda*(Rplanet/rp - 1.0) - 0.5*(v_wp/csp)**2)
            rhoZone = (RHplanet/Pplanet)*presZone
            ekinZone = 0.5*(velxZone*velxZone + velyZone*velyZone + velzZone*velzZone)
            eintZone = presZone/(rhoZone*(sim_gamma - 1.0))
            enerZone = max(eintZone + ekinZone, sim_smallP)
            blkPtr(DENS_VAR,i,j,k) = rhoZone
            blkPtr(PRES_VAR,i,j,k) = presZone
            blkPtr(VELX_VAR,i,j,k) = velxZone
            blkPtr(VELY_VAR,i,j,k) = velyZone
            blkPtr(VELZ_VAR,i,j,k) = velzZone
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
