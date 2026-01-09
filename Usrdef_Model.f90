! Copyright 2013 RBINS-MUMM
!
! Licensed under the EUPL, Version 1.1 or - as soon they will be approved by
! the European Commission - subsequent versions of the EUPL (the "Licence");
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!
! http://ec.europa.eu/idabc/eupl
!
! Unless required by the applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.

!****************************************************************************
!
! *Usrdef_Model* User-defined model setup
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.8
!
! $Date: 2015-10-06 17:42:25 +0200 (Tue, 06 Oct 2015) $
!
! $Revision: 886 $
!
! Description - test case for Red River Coastal area
!
! Reference -May 2016
!
! Subroutines - usrdef_init_params, usrdef_mod_params, usrdef_grid,
!               usrdef_partition, usrdef_phsics, usrdef_1dsur_spec,
!               usrdef_2dobc_spec, usrdef_profobc_spec, usrdef_1dsur_data,
!               usrdef_2dobc_data, usrdef_profobc_data, usrdef_rlxobc_spec
!
!****************************************************************************
!

!============================================================================

SUBROUTINE usrdef_init_params
!****************************************************************************
!
! *usrdef_init_params* Define parameters for monitoring
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description - test case for Red River coastal area
!
! Reference -
!
! Calling program - simulation_start
!
!****************************************************************************
!
USE iopars
USE paralpars

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iproc

!--------------------------
!1. Cold/warm start--------
!--------------------------
cold_start = .false.
! siempre se escribe como false
!--------------------------

!--------------------------
!3. Log files
!--------------------------
!---program leveling in log files----
if (master) then
     levprocs_ini(1:4) = 1  
     levprocs_run(1:4) = 1 ! The size (npworld) of levprocs_run equals either the number of processes defined within the MPI COMM WORLD communicator in the parallel case. Default is zero. 
else
   levprocs_ini = 0 
   levprocs_run = 0
endif

!     levprocs_ini(1:4) = 1  
!levprocs_ini = 1 
!levprocs_run = 1 
!------------------------------------

!-----------------
!6. Timing
!-----------------
levtimer = 3 ! escribe el tiempo para correrlo en modo paralelo. For each part of the simulation (hydrodynamics, 2D mode, etc), “timing” information are given, one for each processor
!-----------------

!--------------------------------
!7------------Parallel setup-----
!--------------------------------
if (npworld.gt.1) nprocscoh = 10 ! si el número de procesos usados por COHERENS por default es mayor que 1, entonces que sea 72. ¿Por qué 72? ! 10 procesadores
!--------------------------------

RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params* Define parameters for physical model
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.8
!
! Description - RRD coastal area
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
USE gridpars
USE iopars
USE paralpars
USE nestgrids
USE physpars
USE switches
USE syspars
USE tide
USE timepars
USE turbpars
USE time_routines, ONLY: convert_date, log_timer_in, log_timer_out

IMPLICIT NONE
INTEGER :: iyear, imont, i, iflag,  iyear_initial, imont_initial 
INTEGER, DIMENSION(7) :: intdate
CHARACTER (LEN=leniofile) :: filename, initialfile


! these lines are to write a final condition each month, so if there is
! any problem and the model crashes, I don't lose everything :(
! final condition once every month ----------------------------
integer, dimension(:), allocatable :: days_months ! Days in months (flexible size)
logical :: is_leap_year               ! Flag for leap year
integer :: syear, l                          ! Loop variable
!  -------------------------------------------------------------

procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()




!-------------------------
!---Initialise parameters-
!-------------------------

!1. Process numbers----------------------
!if (npworld.GT.1) then
!  ---total number of processes
!nprocs = 10
!  ---number of processes in X-direction
!nprocsx = 5
!  ---number of processes in Y-direction
!nprocsy = 2 !se multiplica 2*5
!endif
!----------------------------------------

!--------------------------------------------------------
!-------------------Switches-----------------------------
!--------------------------------------------------------
!---grid dimension (1/2/3)--
iopt_grid_nodim = 3 !Dimension of the model grid. 3-D grid. Default.
iopt_grid_sph = 1 !Type of coordinates. 0: Cartesian coordinates.
iopt_grid_htype = 1 !1 !type of horizontal grid. 1: Uniform rectangular. Default
iopt_grid_vtype = 1
!---------------------------

!---implicit scheme--
iopt_hydro_impl = 1 !1 
!no, selecciono 1 porque..... buscarlo
!--------------------

!---density---------
iopt_sal = 0
iopt_dens = 3 !
sal_ref = 0
iopt_dens_grad = 3 ! 
iopt_dens_convect = 1

!-----Salinity--------
sal_ref = 0 !

!Optical parameters
optattcoef1_cst = 10 !Default
optattcoef2_cst = 2!0.067 !Default
opt_frac = 0.54 ! Default


!----Bottom boundary-----------------
!The 3-D momentum and the depth-integrated momentum equations are integrated using an implicit algorithm which removes the CFL stability constraint on the time step. This avoids a separate integration of the 2-D equations (except for the 2-D continuity equation (5.32) needed for updating the surface elevation).
!---coefficient (0/1/2/3/4)
iopt_bstres_drag = 3 ! Default. Formulation for the botton drag coefficient: Spatially uniform roughness length.

!---bottom stress formulation (0/1/2)
iopt_bstres_form =2 ! Default, recommended. Formulation for the bottom stress. Quadratic bottom stress.
iopt_bstres_nodim = 3 !Type of currents used in the bottom stress formulation. 3D current taken at the bottom grid cell. Default in case of a 3D simulation (recommended)


!----Turbulence------
!iopt_turb_ntrans = 2

!---horizontal diffusion---
iopt_hdif_coef = 2 !Horizontal diffusion is disabled. Default 
iopt_hdif_scal = 1 !Disabled. Default
iopt_hdif_2D = 1 !Default
iopt_hdif_3D = 1 !Default
iopt_hdif_turb = 0


! Turbulent scheme
iopt_vdif_coef = 2!Algebraic scheme. Non Default 
iopt_turb_alg = 2!Tipe of algebraic scheme



!---advection scheme for 2-D currents (0/1/2/3/4)--
iopt_adv_2D = 3!1 es default
!---advection scheme for 3-D currents (0/1/2/3/4)--
iopt_adv_3D = 3!1 Default 
iopt_adv_scal = 3 !Default




!---meteo--------------
iopt_meteo = 1 !non-flux format. I need (Uw , Vw , SP , Ta , RH, TCC , Prec ) variables. Page 159 manual.
iopt_meteo_stres = 1 !enables input of wind
iopt_meteo_data = 1 !non flux format
iopt_meteo_heat = 1 ! enables input of temperature
iopt_obc_invbar = 1
iopt_meteo_precip = 1
iopt_meteo_pres = 1
iopt_sflux_precip = 1
!----------------------



!---temperature-------
!non -default configurations to let temperature change over space and time
iopt_temp = 2!Saheed
iopt_temp_sbc = 1 !Saheed. Default. Neumann surface boundary condition for temperature
iopt_temp_optic = 1 ! Solar radiation is absorbed within the water column using specified values for the attenuation depths. Default (recommended).



!---implicit switches---
iopt_mg_cycle = 1 !Type of cycle used in the multi-grid procedure. 1: V-cycle. Default.
iopt_mg_prolong = 2 !Type of multi-grid prolongation operator. 2: Bilinear interpolation. Default.
iopt_mg_smoother = 1 !Type of smoother for the multi-grid scheme. 1: Jacobi. Default.
!-----------------------

! Surface fluxes
iopt_sflux_qlong  = 2 !Upwards surface long-wave radiation
iopt_sflux_qshort = 3 !Downward surface short-wave radiation

!Estratificacion atmosferica
iopt_sflux_strat = 1 ! Monin-Obukhov formulation
iopt_sflux_pars = 5 ! Kondo atmospheric stratification
!---------------------




!-------------------------

!---wet-dry---------------
!ntobcrlx=2000
!flooding
!iopt_fld = 1
!dcrit_fld = depmin
!-------------------------

!-----netcdf output-----------------------------------------
iopt_CDF_format = 2! sAHEED1!2!64 bit, no file size limit (=1,classic)
iopt_CDF_tlim = 2! SAHEED1!2 ! dimensión temporal ilimitada
!-----------------------------------------------------------

!-------------------------------------------
!3. Date/time parameters
!-------------------------------------------
!read(runtitle(1:6),'(i4,i2)') iyear,imont !The lines above aims to set the year and the month to the first 4 digits and the 5 th & 6 th digit respectively.


!-----------------------------------------------
CStartDateTime(1:19) = '2004/01/01;00:00:00'!---
CEndDateTime(1:19) = '2004/06/01;00:00:00'  !---
!-----------------------------------------------

!-------------------time step----------------------------------
timestep =60.0 !MERGE(2.0,60.0,iopt_hydro_impl.eq.0)
!---counter for 3-D mode
ic3d = 1! 10!1!MERGE(1,900,iopt_hydro_impl.eq.1.OR.iopt_grid_nodim.eq.2)
!--------------------------------------------------------------

!------------Physical model constants--------------------------

!--------------grid dimensions---------
nc = 145 +1 ; nr = 105 +1
nz = 10
!--------------------------------------

!---number of open sea boundaries------
nosbv = 0
nosbu = 0
!--------------------------------------

!---number of open river boundaries----
!nrvbu = 0; nrvbv = 0
!--------------------------------------


!---uniform bottom roughness length----
zrough_cst =0.5! 0.25!0.0035 !
specheat = 4184
!--------------------------------------

!---multigrid parameters---------------
nomglevels = 1!4!MERGE(4,3,iopt_MPI.eq.0)
mg_tol = 1.0E-07 !Default COHERENS          !1.0E-04 Saheed
ur_smooth = 0.8 ! este es el default y el usado por Saheed
!------------------------------------------------------------------

!------------------------tidal indices---------------------------------------------------- 
!index_obc(1:nconobc) = (/icon_M2, icon_S2, icon_N2, icon_K2, icon_K1, icon_O1, icon_P1, & 
!                        & icon_Q1, icon_MF, icon_MM, icon_M4, icon_MS4, icon_MN4/)
!-----------------------------------------------------------------------------------------

!--------------Model I/O file properties--------------------------------

!----------------------- Input------------------------------------------
!----model grid---------------------------------------------------------
modfiles(io_modgrd,1,1)%status = 'N'!User-defined in which case a corresponding usrdef routine is called.
modfiles(io_modgrd,1,1)%form = 'A'
modfiles(io_modgrd,1,1)%filename ='DATA/Lanalhue_grid3.dat'!Lanalhue_grid.dat' 
!-----------------------------------------------------------------------

!---------open boundary conditions (2-D)--------------------------------
!modfiles(io_2uvobc,1,1)%status = '0'
!modfiles(io_2uvobc,1,1)%form = 'A'
!modfiles(io_2uvobc,1,1)%filename = 'DATA/hlb_coastal_obc.dat'
!-----------------------------------------------------------------------
!modfiles(io_2uvobc,2,1)%status = 'N'
!modfiles(io_2uvobc,2,1)%tlims = (/0,0,1/)

!---restart times (needed for creating initial conditions)--------------
norestarts = 1
ntrestart(1) = int_fill 
initialfile = outtitle
!if (imont.ne.1) then
!   imont_initial= imont -1
!   write(initialfile(1:6),'(i4.4,i2.2)') iyear,imont_initial
!else
!   imont_initial = 12
!   iyear_initial = iyear - 1
!   write(initialfile(1:6),'(i4.4,i2.2)') iyear_initial,imont_initial
!endif
!------------------------------------------------------------------------

!------initial conditions------------------------------------------------
 ! Temperature profile  
   !modfiles(io_inicon,1,1)%status = 'N'
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   !modfiles(io_inicon,1,1)%form = 'A'
   modfiles(io_inicon,ics_phys,1)%form = 'N'
   !modfiles(io_inicon,1,1)%filename = 'DATA/initcon/ll_HD_init_112001_10lev.dat'
   modfiles(io_inicon,ics_phys,1)%filename ='/home/maria/my_models4/DATA/initcon/ll_init_jan2004_physics.nc'

! velocidades por defecto iguales a cero

! Salinidad, temperatura y turbulencia. Por ahora no necesito

modfiles(io_fincon,1,2)%status = 'W'
modfiles(io_fincon,1,2)%form = 'N'
modfiles(io_fincon,1,2)%filename = 'DATA/initcon/best_bugHugo/ll_init_jun2004_physics.nc'
!........................................................................















!-------------------meteo------------------------------------------------
modfiles(io_metsur,1,1)%status = 'N'
modfiles(io_metsur,1,1)%filename = '/home/maria/my_models4/DATA/meteo/era5/era5_2000_2010.nc'
!modfiles(io_metsur,1,1)%filename = 'DATA/era5_ll_2001_2003_rate.nc'
modfiles(io_metsur,1,1)%form = 'N'
modfiles(io_metsur,1,1)%tlims = (/0,int_fill,60/)
!------------------------------------------------------------------------

!---open boundary conditions (3-D current)-------------------------------
!modfiles(io_3uvobc,1:2,1)%status = 'N'
!modfiles(io_3uvobc,2,1)%tlims = (/0,0,1/)
!------------------------------------------------------------------------

!--------open boundary conditions (salinity)-----------------------------
!------------------------------------------------------------------------

!7. Surface grid parameters
!--------------------------
surfacegrids(igrd_model,1)%x0dat = -73.39 -0.0005 !0.00025!0.0005
surfacegrids(igrd_model,1)%y0dat = -37.99 -0.0005 !0.00025!0.0005
surfacegrids(igrd_model,1)%delxdat = 0.001 !0.0005! 0.001
surfacegrids(igrd_model,1)%delydat = 0.001 !0.0005! 0.001

!--------meteo grid parameters-------------------------------------------
!--------meteo grid parameters-------------------------------------------
!esto debe estar activo si quiero correr el modelo con datos meteorologicos
!ellos son el unico forzante en mi lago, puesto que no tengo oceano`
surfacegrids(igrd_meteo,1)%nhtype = 1 !(grid descriptor, file number). file numer es 1 si se esta usando para input
surfacegrids(igrd_meteo,1)%n1dat = 13 !nc of wind data
surfacegrids(igrd_meteo,1)%n2dat = 13 !nr of wind data
surfacegrids(igrd_meteo,1)%x0dat = -74.5-0.125  
surfacegrids(igrd_meteo,1)%y0dat = -39.5-0.125
surfacegrids(igrd_meteo,1)%delxdat = 0.25
surfacegrids(igrd_meteo,1)%delydat = 0.25
!-----------------------time series output-------------------------------
nosetstsr = 2
novarstsr = 6
!------------------------------------------------------------------------

CALL log_timer_out()
RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define model grid arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - test case bohai
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - close_filepars, open_filepars
!
!************************************************************************
USE depths
USE grid
USE gridpars
USE iopars
USE switches
USE inout_routines, ONLY: close_filepars, open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ii, i, iunit, j, jj
REAL :: depmin, depmean_flag


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!1. Open grid file
!-----------------
!

CALL open_filepars(modfiles(io_modgrd,1,1))
iunit = modfiles(io_modgrd,1,1)%iunit

!
!2. Water depths
!---------------



!depmean_flag = -999.9999

READ (iunit,*)
READ (iunit,*)
j_210: DO j=1,nr-1
   READ (iunit,*) depmeanglb(1:nc-1,j)
ENDDO j_210

depmin = 0.1!18 ! esto al ser menor que 2, si hubieran mareas, tiraria error. Pero cocmo aca es un lago y no hay mareas, da lo mismo.
j_220: DO j=1,nr-1
i_220: DO i=1,nc-1
       IF (depmeanglb(i,j).le.0.0) THEN
          depmeanglb(i,j) = 0
       ELSE
          depmeanglb(i,j) = MAX(depmin,depmeanglb(i,j))
       ENDIF
ENDDO i_220
ENDDO j_220




!!3. Open boundary locations

!--------------------------
!
!---U-nodes

CALL close_filepars(modfiles(io_modgrd,1,1))

CALL log_timer_out()
!--------------------------------------------------------------
RETURN
END SUBROUTINE usrdef_grid

!========================================================================

SUBROUTINE usrdef_partition
!************************************************************************
!
! *usrdef_partition* Define domain decompositionnd
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - domain_decomposition
!
!************************************************************************

IMPLICIT NONE

RETURN

END SUBROUTINE usrdef_partition

!========================================================================

SUBROUTINE usrdef_phsics
!************************************************************************
!
! *usrdef_phsics* Define initial conditions for physics
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE density
USE gridpars
USE iopars
USE error_routines, ONLY: nerrs, error_abort
USE inout_routines, ONLY: close_filepars, open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


!
!----------Local variables------------------------------
INTEGER :: iunit, k, nnz
REAL, DIMENSION(nz) :: tin
!-------------------------------------------------------
procname(pglev+1) = 'usrdef_phsics'
CALL log_timer_in()
!-------------------------------------------------------
!---------Temperature-----------------------------------
!---open file with initial profile
! esto es solo para leer el perfil inicial
! solo cuando estoy leyendo un archivo en un formato externo a coherens

!CALL open_filepars(modfiles(io_inicon,1,1))
!iunit = modfiles(io_inicon,1,1)%iunit

!-----------read initial conditions---------------------
!read (iunit,*) tin
!k_310: do k=1,nz
!   temp(1:ncloc,1:nrloc,k) = tin(k)
!enddo k_310
!CALL close_filepars(modfiles(io_inicon,1,1))


CALL log_timer_out()

RETURN
END SUBROUTINE usrdef_phsics

!========================================================================

SUBROUTINE usrdef_1dsur_spec
!************************************************************************
!
! *usrdef_1dsur_spec* Define specifier arrays for 1-D forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_1dsur_spec
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_1dsur_spec

!========================================================================

SUBROUTINE usrdef_2dobc_spec(iddesc)
!************************************************************************
!
! *usrdef_2dobc_spec* Define specifier arrays for 2-D open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case bohai
!
! Reference -
!
! Calling program - define_2dobc_spec
!
! Module calls - close_file, open_file 
!
!************************************************************************
!
USE iopars
USE obconds
USE syspars
USE gridpars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE inout_routines, ONLY: close_filepars, open_filepars, close_file, open_file
IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc
!*Local variables
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iunit, kk, i, j, l
REAL :: dummy1, dummy2, dummy3, dummy4
PRINT *, 'In 2dobc_spec'
procname(pglev+1) = 'usrdef_2dobc_spec'

CALL log_timer_in()

CALL log_timer_out()
RETURN
!-----------------------------------------------------------------------
END SUBROUTINE usrdef_2dobc_spec

!=======================================================----=================

SUBROUTINE usrdef_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                             & noprofsd,indexprof,nofiles,nobux,nobvy,novars)
!****************************************************************************
!
! *usrdef_profobc_spec* Define specifier arrays for open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.x
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_profobc_spec
!
!*****************************************************************************
!
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, nobux, nobvy, nofiles, novars
INTEGER, INTENT(INOUT), DIMENSION(2:nofiles,novars) :: noprofsd
INTEGER, INTENT(INOUT), DIMENSION(nobux) :: itypobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy) :: itypobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux,novars) :: iprofobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy,novars) :: iprofobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux+nobvy,2:nofiles,novars) :: indexprof
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id
!*itypobux*  INTEGER Type of U- or X-open boundary condition
!*itypobvy*  INTEGER Type of V- or Y-open boundary condition
!*iprofobux* INTEGER Profile numbers at U- or X-open boundaries
!*iprofobvy* INTEGER Profile numbers at V- or Y-open boundaries
!*noprofsd*  INTEGER Number of profiles per data file
!*indexprof* INTEGER Mapping array (for each data variable) of the profile
!                    numbers in the data files to the profile numbers assigned
!                    to the open boundaries. The physical size of the first
!                    dimension equals the number of profiles in a data file.
!*nofiles*   INTEGER Number of data files (+1)
!*nobux*     INTEGER Number of nodes at U- or X-open boundaries
!*nobvy*     INTEGER Number of nodes at V- or Y-open boundaries
!*novars*    INTEGER Number of data variables
!
!------------------------------------------------------------------------------
!local variables
INTEGER :: kk

procname(pglev+1) = 'usrdef_profobc_spec'
CALL log_timer_in()


CALL log_timer_out()

RETURN

END SUBROUTINE usrdef_profobc_spec

!========================================================================
SUBROUTINE usrdef_1dsur_data(ciodatetime,data1d,novars)
!************************************************************************
!
! *usrdef_1dsur_data* Define data for 1-D forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_1dsur_data
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: novars
REAL, INTENT(INOUT), DIMENSION(novars) :: data1d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time in data file
!*data1d*      REAL    Surface data
!*novars*      INTEGER Number of surface data
!
!------------------------------------------------------------------------------
RETURN
END SUBROUTINE usrdef_1dsur_data
!========================================================================
SUBROUTINE usrdef_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *usrdef_2dobc_data* Define open boundary data for 2-D mode
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description -
!
! Reference -
!
! Calling program - define_2dobc_data
!
!************************************************************************
!
USE gridpars
USE iopars
USE syspars
USE grid
USE timepars
USE time_routines, ONLY: add_secs_to_date, convert_date, log_timer_in, log_timer_out
USE inout_routines, ONLY: close_file, open_file, open_filepars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: data2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*data2d*      REAL    Input data
!*nodat*       INTEGER Number of input data
!*novars*      INTEGER Number of input parameters
!
!------------------------------------------------------------------------------
procname(pglev+1) = 'usrdef_2dobc_data'
CALL log_timer_in()

CALL log_timer_out()
RETURN
END SUBROUTINE usrdef_2dobc_data

!========================================================================

SUBROUTINE usrdef_profobc_data_3d(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
!************************************************************************
!
! *usrdef_profobc_data* Define physical open boundary profiles
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.2
!
! Description -
!
! Reference -
!
! Calling program - define_profobc_data
!
!************************************************************************
!
USE gridpars
USE iopars  
USE syspars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out, convert_date,add_secs_to_date
USE inout_routines, ONLY: close_file, open_file


IMPLICIT NONE

!*  Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, numprofs
REAL, INTENT(INOUT), DIMENSION(numprofs,nz) :: psiprofdat

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*psiprofdat*  REAL     Profile arrays
!*numprofs*    INTEGER  Number of profiles in data file
!*nobcvars*    INTEGER  Effective number of data variables for which open
!                       boundary conditions are applied
!
!------------------------------------------------------------------------------
!
! local variables
procname(pglev+1) = 'usrdef_profobc_data_3d'
CALL log_timer_in()

CALL log_timer_out()

RETURN

END SUBROUTINE usrdef_profobc_data_3d

!========================================================================
