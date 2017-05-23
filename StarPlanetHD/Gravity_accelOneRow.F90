!!****if* source/physics/Gravity/GravityMain/PointMass/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow 
!!
!! SYNOPSIS
!!
!!  call Gravity_accelOneRow(integer(IN)  :: pos(2),
!!                           integer(IN)  :: sweepDir,
!!                           integer(IN)  :: blockID,
!!                           integer(IN)  :: numCells,
!!                           real(INOUT)  :: grav(numCells),
!!                           integer(IN),optional :: potentialIndex,
!!                           integer(IN),optional :: extraAccelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration for a row
!!  of cells in a specified direction in a given block.
!!
!! ARGUMENTS
!!
!!  pos      :  Row indices transverse to the sweep direction
!!  sweepDir :    The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y, and SWEEP_Z. These values are defined
!!              in constants.h.
!!  blockID  :  The local identifier of the block to work on
!!  numCells :  Number of cells to update in grav()
!!  grav()   :   Array to receive result
!!  potentialIndex :  optional, not applicable in pointmass gravity
!!  extraAccelVars :  optional, ignored in this implementation
!! 
!!***

subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, &
                                potentialIndex, extraAccelVars)

!=======================================================================

  use Simulation_data
  use Gravity_data, ONLY: grv_ptxpos, grv_ptypos, grv_ptzpos, grv_factor, &
       useGravity
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getBlkPtr

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: sweepDir,blockID,numCells
  integer, dimension(2),INTENT(in) ::pos
  real, dimension(numCells),INTENT(inout) :: grav
  integer,intent(IN),optional :: potentialIndex
  integer,intent(IN),OPTIONAL :: extraAccelVars(MDIM)

!==========================================================================


#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_KHI_GC) :: zCenter
  real,dimension(GRID_JHI_GC) :: yCenter
  real,dimension(GRID_IHI_GC) :: xCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
  integer, dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
#endif
  real :: dr32, tmpdr32

  integer :: sizeX,sizeY,sizez

  integer :: ii,j,k
  logical :: gcell = .true.

  real, pointer, dimension(:,:,:,:) :: solnData 

!==============================================================================

  real :: Big_G, Rsun, Msun, Mstar, Rstar, RHstar, pi
  real :: Rj, Mj, Mplanet, Rplanet, RHplanet, rp, rp_c
  real :: gfac, gfac_int, gfacp, gfac_intp, a, omega_orb, omegas, omegap, omega_fr
  real :: Fin_x, Fin_y, vx, vy
  ! Simulation units.
  real :: UNIT_DENSITY, UNIT_LENGTH, UNIT_VELOCITY
  real :: UNIT_MP, UNIT_G, UNIT_TIME, UNIT_MASS

  pi = 3.1415926536
  Big_G = 6.67259e-8
  Msun = 1.989e+33    ![g]
  Rsun = 6.955e+10    ![cm]
  Mj = 0.0009543*Msun ![g]
  Rj = 0.10045*Rsun   ![cm]

  UNIT_DENSITY = 1.0e-15
  UNIT_LENGTH = 6.955e+10
  UNIT_VELOCITY = 1.0e+5
  UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY
  UNIT_MASS = UNIT_DENSITY*UNIT_LENGTH**3
  UNIT_G = Big_G*UNIT_DENSITY*UNIT_TIME**2

  Rstar = sim_Rs*Rsun/UNIT_LENGTH
  Mstar = sim_Ms*Msun/UNIT_MASS
  RHstar = sim_RHs/UNIT_DENSITY

  Rplanet = sim_Rp*Rj/UNIT_LENGTH
  Mplanet = sim_Mp*Mj/UNIT_MASS
  RHplanet = sim_RHp/UNIT_DENSITY

  a = sim_sep*0.01*1.496e+13/UNIT_LENGTH
  omega_orb = sqrt(UNIT_G*Mstar/(a**3))
  omegas    = omega_orb
  omegap    = omega_orb
  omega_fr  = omega_orb

  gfac = -UNIT_G*Mstar
  gfac_int = -(4.0/3.0)*pi*UNIT_G*RHstar

  gfacp = -UNIT_G*Mplanet
  gfac_intp = -(4.0/3.0)*pi*UNIT_G*RHplanet

!==============================================================================
  if (.NOT.useGravity) return

  j=pos(1)
  k=pos(2)
#ifndef FIXEDBLOCKSIZE
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)
  sizeY=blkLimitsGC(HIGH,JAXIS)
  sizeZ=blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX))
  allocate(yCenter(sizeY))
  allocate(zCenter(sizeZ))
#else
  sizeX=GRID_IHI_GC
  sizeY=GRID_JHI_GC
  sizeZ=GRID_KHI_GC
#endif
  zCenter = 0.
  yCenter = 0.
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
     !zCenter = zCenter - grv_ptzpos
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCenter, sizeY)
     !yCenter = yCenter - grv_ptypos
  endif
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)
  !xCenter = xCenter - grv_ptxpos
 
  call Grid_getBlkPtr(blockID,solnData,CENTER) 

  if (sweepDir .eq. SWEEP_X) then                       ! x-component

     tmpdr32 = yCenter(j)*yCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells

        dr32 = sqrt(xCenter(ii)*xCenter(ii) + tmpdr32)
        dr32 = dr32*dr32*dr32

        rp = sqrt((xCenter(ii)-a)*(xCenter(ii)-a) + tmpdr32)
        rp = rp*rp*rp

        if (sqrt(xCenter(ii)*xCenter(ii) + tmpdr32) > Rstar) then
          vy = solnData(VELY_VAR,ii,j,k)
          Fin_x = omega_fr*omega_fr*xCenter(ii) + 2.0*omega_fr*vy
          grav(ii) = gfac*xCenter(ii)/dr32 + gfacp*(xCenter(ii)-a)/rp + Fin_x
        else if (sqrt(xCenter(ii)*xCenter(ii) + tmpdr32) < Rstar) then
          grav(ii) = gfac_int*xCenter(ii) + gfac_intp*(xCenter(ii)-a)
        end if

     enddo


  else if (sweepDir .eq. SWEEP_Y) then          ! y-component

     tmpdr32 = xCenter(j)*xCenter(j) + zCenter(k)*zCenter(k) 

     rp_c = (xCenter(j)-a)*(xCenter(j)-a) + zCenter(k)*zCenter(k)

     do ii = 1, numCells
        
        dr32 = sqrt(yCenter(ii)*yCenter(ii) + tmpdr32)
        dr32 = dr32*dr32*dr32

        rp = sqrt(yCenter(ii)*yCenter(ii) + tmpdr32)
        rp = rp*rp*rp

        if (sqrt(yCenter(ii)*yCenter(ii) + tmpdr32) > Rstar) then
          vx = solnData(VELX_VAR,j,ii,k)
          Fin_y = omega_fr*omega_fr*yCenter(ii) - 2.0*omega_fr*vx
          grav(ii) = gfac*yCenter(ii)/dr32 + gfacp*(yCenter(ii))/rp + Fin_y
        else if (sqrt(yCenter(ii)*yCenter(ii) + tmpdr32) < Rstar) then
          grav(ii) = gfac_int*yCenter(ii) + gfac_intp*yCenter(ii)
        end if

     enddo

  else if (sweepDir .eq. SWEEP_Z) then          ! z-component

     tmpdr32 = xCenter(j)*xCenter(j) + yCenter(k)*yCenter(k) 

     rp_c = (xCenter(j)-a)*(xCenter(j)-a) + yCenter(k)*yCenter(k)
     rp = rp*rp*rp

     do ii = 1, numCells
        
        dr32 = sqrt(zCenter(ii)*zCenter(ii) + tmpdr32)           
        dr32 = dr32*dr32*dr32

        rp = sqrt(zCenter(ii)*zCenter(ii) + tmpdr32)
        
        if (sqrt(zCenter(ii)*zCenter(ii) + tmpdr32) > Rstar) then
          grav(ii) = gfac*zCenter(ii)/dr32 + gfacp*(zCenter(ii))/rp
        else if (sqrt(zCenter(ii)*zCenter(ii) + tmpdr32) < Rstar) then
          grav(ii) = gfac_int*zCenter(ii) + gfac_int*zCenter(ii)
        end if

     enddo

  endif

!==============================================================================
#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
#endif

  return

end subroutine Gravity_accelOneRow
