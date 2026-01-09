!
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

!************************************************************************
!
! *Usrdef_Surface_Data* User-defined surface data setup
!
! Author - Romanelli Hugo
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V3.0
!
! Date: 2024-06-01
!
! Description - Belgian Coast (BeC)
!
! Subroutines - usrdef_surface_absgrd, usrdef_surface_relgrd,
!               usrdef_surface_data
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_surface_absgrd(iddesc,ifil,n1dat,n2dat,xcoord,ycoord)
!************************************************************************
!
! *usrdef_surface_absgrd* Define coordinate arrays of surface grid(s)
!                         with respect to model grid
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
!************************************************************************
!
USE datatypes
USE gridpars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, n1dat, n2dat
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat) :: xcoord, ycoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Grid file id
!*ifil*      INTEGER No. of grid file
!*n1dat*     INTEGER X-dimension of data grid
!*n2dat*     INTEGER Y-dimension of data grid
!*xcoord*    REAL    X-coordinates of data grid
!*ycoord*    REAL    Y-coordinates of data grid
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_surface_absgrd

!========================================================================

SUBROUTINE usrdef_surface_relgrd(iddesc,ifil,surfgridglb,nx,ny)
!************************************************************************
!
! *usrdef_surface_relgrd* Define relative coordinate array of surface grid(s)
!                         with respect to model grid
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
!************************************************************************
!
USE datatypes
USE gridpars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, nx, ny
TYPE (HRelativeCoords), INTENT(INOUT), DIMENSION(nx,ny) :: surfgridglb

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER Grid file id
!*ifil*        INTEGER No. of grid file
!*surfgridglb* DERIVED Relative coordinates of model grid to data grid or data
!                      grid to model grid
!*nx*          INTEGER X-dimension of data grid
!*ny*          INTEGER Y-dimension of data grid 
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_surface_relgrd

!========================================================================

SUBROUTINE usrdef_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,&
                             & novars)
!************************************************************************
!
! *usrdef_surface_data* Define surface input data
!
! Calling program - define_surface_data
!
!************************************************************************
!
USE iopars
USE netcdf
USE syspars
USE time_routines

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, novars, n1dat, n2dat
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat,novars) :: surdata

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!                io_metsur => meteo data
!                io_sstsur => sea surface temperature data
!                io_wavsur => wave data
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*surdata*     REAL     Data array
!*n1dat*       INTEGER  X-dimension of data array
!*n2dat*       INTEGER  Y-dimension of data array
!*novars*      INTEGER  Number of data parameters
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL, SAVE :: FIRST = .TRUE.
INTEGER :: ncid, status
INTEGER :: ntimeid, timeid
INTEGER, SAVE :: ntime
INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: time
INTEGER, SAVE :: u10id, v10id, mslid,t2mid, tccid, d2mid, avg_tprateid
INTEGER, DIMENSION(7), SAVE :: refdate
INTEGER, SAVE :: record = 0
INTEGER, DIMENSION(7) :: intdate
REAL, DIMENSION(n1dat, n2dat) :: u10, v10, msl, t2m, tcc, d2m, d2m_tot, avg_tprate
REAL, DIMENSION(n1dat, n2dat) :: numerator, denominator 
REAL :: alpha, beta


!---Initialization on first call

IF (FIRST) THEN

  ! Open file
  status = nf90_open(modfiles(io_metsur, 1, 1)%filename, nf90_nowrite, ncid)
  modfiles(io_metsur, 1, 1)%iostat = 1

  ! Time
  status = nf90_inq_dimid(ncid, 'valid_time', ntimeid)
  status = nf90_inquire_dimension(ncid, ntimeid, len = ntime)
  status = nf90_inq_varid(ncid, 'valid_time', timeid)
  allocate(time(ntime))
  status = nf90_get_var(ncid, timeid, time)

  ! Meteo variables
  status = nf90_inq_varid(ncid, 'u10', u10id)
  status = nf90_inq_varid(ncid, 'v10', v10id)
  status = nf90_inq_varid(ncid, 'msl', mslid)
  status = nf90_inq_varid(ncid, 't2m', t2mid)
  status = nf90_inq_varid(ncid, 'd2m', d2mid)
  status = nf90_inq_varid(ncid, 'tcc', tccid)
  status = nf90_inq_varid(ncid, 'avg_tprate', avg_tprateid)

  ! Close file
  status = nf90_close(ncid) 

  ! Reference date
  refdate = (/1970, 1, 1, 0, 0, 0, 0/)

  ! End of initialization
  FIRST = .FALSE.


!---Read/update data
ELSE

  ! Update record
  record = record + 1
  IF (record > ntime) THEN
    write(*, *) 'End of meteo file reached'
  ENDIF

  ! Record time
  CALL add_secs_to_date_int(refdate, intdate, time(record), 1.)
  ciodatetime = convert_date(intdate)

  ! Open file
  status = nf90_open(modfiles(io_metsur, 1, 1)%filename, nf90_nowrite, ncid)

  ! Read data
  status = nf90_get_var(ncid, u10id, u10, start = (/1, 1, record/), count = (/n1dat, n2dat, 1/))
  status = nf90_get_var(ncid, v10id, v10, start = (/1, 1, record/), count = (/n1dat, n2dat, 1/))
  status = nf90_get_var(ncid, mslid, msl, start = (/1, 1, record/), count = (/n1dat, n2dat, 1/))
  status = nf90_get_var(ncid, t2mid, t2m, start = (/1, 1, record/), count = (/n1dat, n2dat, 1/))
  status = nf90_get_var(ncid, d2mid, d2m, start = (/1, 1, record/), count = (/n1dat, n2dat, 1/))
  status = nf90_get_var(ncid, tccid, tcc, start = (/1, 1, record/), count = (/n1dat, n2dat, 1/))
  status = nf90_get_var(ncid, avg_tprateid, avg_tprate, start = (/1, 1, record/), count = (/n1dat, n2dat, 1/))

  ! Close file
  status = nf90_close(ncid)

  ! Store data (u10, v10, msl, temp, RH, fc, avg_tprate)
  surdata(:, :, 1) = u10(:, n2dat:1:-1)
  surdata(:, :, 2) = v10(:, n2dat:1:-1)
  surdata(:, :, 3) = msl(:, n2dat:1:-1)
  surdata(:, :, 4) = t2m(:, n2dat:1:-1) - 273.15
  
  !---begin computing relative humidity with Magnus formula
  alpha = 243.04 !°C
  beta = 17.625
  d2m_tot = d2m(:, n2dat:1:-1) - 273.15
  numerator = exp((beta*d2m_tot)/(alpha+d2m_tot))
  denominator = exp((beta*surdata(:,:,4))/(alpha+surdata(:,:,4)))
  surdata(:,:,5) = (numerator/denominator)
  !---end computing relative humidity

  surdata(:,:,6) = tcc(:, n2dat:1:-1)
  surdata(:,:,7) = avg_tprate(:, n2dat:1:-1)
  
ENDIF

RETURN

END SUBROUTINE usrdef_surface_data
