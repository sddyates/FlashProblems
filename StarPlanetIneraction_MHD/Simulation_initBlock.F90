!!****if* source/Simulation/SimulationMain/magnetoHD/Rotor/Simulation_initBlock
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
!!  Reference: Balsara and Spicer, JCP, 149:270--292, 1999
!!             Balsara, The Astrophys Suppl Series, 151:148--184, 2004
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  
!!
!! 
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
  real :: enerZone, ekinZone, eintZone
  real, allocatable, dimension(:) :: xCoord,yCoord,zCoord
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData

  ! Constants.
  real :: PI, Big_G, Msun, Rsun, Mj, Rj
  
  ! Geometric quantities.
  real :: x, y, z, rs, thetas, phis

  ! Stellar parameters.
  real :: Ts, B0s, Ms, Rs, RHs

  ! Stellar wind parameters.
  real :: v_esc, css, rcs, v_ws, v_ws_x, v_ws_y, v_ws_z
  real :: RHOs, PRs, PRSs 

  ! Variables used in the Newton-Raphson routine.
  real :: lambda, eta, b, psi, step, f, df
  integer :: o

  ! Simulation units.
  real :: UNIT_DNESIY, UNIT_LENGTH, UNIT_VELOCITY
  real :: UNIT_G, UNIT_TIME, UNIT_MASS, UNIT_MAG

  PI = 3.1415926536
  Big_G = 6.67259e-8  ![cm3 g-1 s-2]
  Msun = 1.989e+33    ![g]
  Rsun = 6.955e+10    ![cm]
  Mj = 0.0009543*Msun ![g]
  Rj = 0.10045*Rsun   ![cm]

  UNIT_TIME = UNIT_LENGTH/UNIT_VELOCITY
  UNIT_MASS = UNIT_DENSITY*UNIT_LENGTH**3
  UNIT_MAG = sqrt(4.0*PI*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)
  UNIT_G = Big_G*UNIT_DENSITY*UNIT_TIME**2

  Ts = sim_Ts
  B0s = sim_B0s/UNIT_MAG
  Ms = sim_Ms*Msun/UNIT_MASS
  Rs = sim_Rs*Rsun/UNIT_LENGTH
  RHs = sim_RHs/UNIT_DENSITY

  v_esc = sqrt(2.0*UNIT_G*Ms/Rs)
  css = sqrt((2.0*Ts)/KELVIN)
  rcs = UNIT_G*Ms/(2.0*css*css)

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

  call Grid_getBlkPtr(blockID, solnData, CENTER)

  #if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID, facexData, FACEX)
     call Grid_getBlkPtr(blockID, faceyData, FACEY)
     call Grid_getBlkPtr(blockID, facezData, FACEZ)
  endif
  #endif

  ! Loop over cells in the block.
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
    z = zCoord(z)
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

        if (r < 0.5*Rs) then

          solnData(DENS_VAR,i,j,k) = Ps + (2.0/3.0)*PI*UNIT_G*RHs*RHs*(Rs*Rs - rs*rs)
          solnData(PRES_VAR,i,j,k) = RHs

          solnData(VELX_VAR,i,j,k) = 0.0
          solnData(VELY_VAR,i,j,k) = 0.0
          solnData(VELZ_VAR,i,j,k) = 0.0

          solnData(MAGX_VAR,i,j,k) = 0.0
          solnData(MAGY_VAR,i,j,k) = 0.0
          solnData(MAGZ_VAR,i,j,k) = 16.0*B0s
          #if NFACE_VARS > 0
          if (sim_killdivb) then
            facexData(MAG_FACE_VAR,i,j,k) = 0.0
            faceyData(MAG_FACE_VAR,i,j,k) = 0.0
            facezData(MAG_FACE_VAR,i,j,k) = 16.0*B0s
          endif
          #endif


        endif

        if (r >= 0.5*Rs .and. r < Rs) then

          solnData(DENS_VAR,i,j,k) = Ps + (2.0/3.0)*PI*UNIT_G*RHs*RHs*(Rs*Rs - rs*rs)
          solnData(PRES_VAR,i,j,k) = RHs

          solnData(VELX_VAR,i,j,k) = 0.0
          solnData(VELY_VAR,i,j,k) = 0.0
          solnData(VELZ_VAR,i,j,k) = 0.0

          solnData(MAGX_VAR,i,j,k) = 3.0*x*z*B0s*Rs**3*rs**-5 
          solnData(MAGY_VAR,i,j,k) = 3.0*z*y*B0s*Rs**3*rs**-5 
          solnData(MAGZ_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0s*Rs**3*rs**-5 

          #if NFACE_VARS > 0
          if (sim_killdivb) then
            facexData(MAG_FACE_VAR,i,j,k) = 3.0*x*z*B0s*Rs**3*rs**-5 
            faceyData(MAG_FACE_VAR,i,j,k) = 3.0*z*y*B0s*Rs**3*rs**-5  
            facezData(MAG_FACE_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0s*(Rs**3*rs**-5  
          endif
          #endif

        endif

        ! The wind external to the star
        if (r >= Rs) then


          ! Newton raphson method for the stellar wind.
          ! This code solves Parkers equation for an 
          ! isothermal stellar wind (Parker 1958).
          lambda = 0.5*(v_escs/css)**2
          eta = rs/Rs
          b = -3.0 - 4.0*log(lambda/2.0) + 4.0*log(eta)+2.0*(lambda/eta)
          o = 0

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
          if isnan(v_ws) then
            v_ws = 0.0
          endif

          v_ws_x = sin(thetas)*(v_ws*cos(phis) + sin(phis)*rs*(omega_fr + omegas))
          v_ws_y = sin(thetas)*(v_ws*sin(phis) - cos(phis)*rs*(omega_fr + omegas))
          v_ws_z = v_ws*cos(thetas)
          PRSs = Ps*exp(lambda*(Rs/rs - 1.0) - 0.5*(v_ws/css)**2)
          RHOs = (RHs/Ps)*PRSs

          solnData(DENS_VAR,i,j,k) = RHOs
          solnData(PRES_VAR,i,j,k) = PRSs

          solnData(VELX_VAR,i,j,k) = v_wp_x
          solnData(VELY_VAR,i,j,k) = v_wp_y
          solnData(VELZ_VAR,i,j,k) = v_wp_z

          solnData(MAGX_VAR,i,j,k) = 3.0*x*z*B0s*Rs**3*rs**-5 
          solnData(MAGY_VAR,i,j,k) = 3.0*z*y*B0s*Rs**3*rs**-5 
          solnData(MAGZ_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0s*Rs**3*rs**-5 

          #if NFACE_VARS > 0
          if (sim_killdivb) then
            facexData(MAG_FACE_VAR,i,j,k) = 3.0*x*z*B0s*Rs**3*rs**-5 
            faceyData(MAG_FACE_VAR,i,j,k) = 3.0*z*y*B0s*Rs**3*rs**-5  
            facezData(MAG_FACE_VAR,i,j,k) = (3.0*z*z - rs*rs)*B0s*(Rs**3*rs**-5  
          endif
          #endif

        endif

        magprs = 0.5* dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                    solnData(MAGX_VAR:MAGZ_VAR,i,j,k)) 
        ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                       solnData(VELX_VAR:VELZ_VAR,i,j,k))
        eintZone = solnData(PRES_VAR,i,j,k)/(sim_gamma - 1.0) & 
                                             /solnData(DENS_VAR,i,j,k)
        enerZone = max(eintZone + ekinZone, sim_smallP)

        solnData(DIVB_VAR,i,j,k) = 0.0
        solnData(MAGP_VAR,i,j,k) = magprs
        solnData(ENER_VAR,i,j,k) = enerZone 
        solnData(EINT_VAR,i,j,k) = eintZone
        solnData(GAMC_VAR,i,j,k) = sim_gamma
        solnData(GAME_VAR,i,j,k) = sim_gamma

      enddo
    enddo
  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID, solnData, CENTER)

  #if NFACE_VARS > 0
  !if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID, facexData, FACEX)
     call Grid_releaseBlkPtr(blockID, faceyData, FACEY)
     call Grid_releaseBlkPtr(blockID, facezData, FACEz)
  !endif
  #endif

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock
