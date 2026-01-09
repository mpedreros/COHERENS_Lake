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
! *Usrdef_Time_Series* Define time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.8
!
! $Date: 2015-05-28 12:01:13 +0200 (Thu, 28 May 2015) $
!
! $Revision: 867 $
!
! Description - test case bohai
!
! Reference -
!
! Routines - usrdef_tsr_params, usrdef_tsr0d_vals, usrdef_tsr2d_vals,
!            usrdef_tsr3d_vals
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_tsr_params
!************************************************************************
!
! *usrdef_tsr_params* Specifiers for time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.8
!
! Description - test case bohai
!
! Reference -
!
! Calling program - time_series_init
!
! Module calls - inquire_var
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE sedids
USE switches
USE timepars
USE modvars_routines, ONLY: inquire_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTEGER :: timeforcheck
INTEGER :: imont, iyear
timeforcheck = nstep
PRINT *, nstep
procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------

!tsrvars%ivarid = (/iarr_uwindatc,iarr_vwindatc,iarr_zeta,&
!                 &iarr_uvel,iarr_vvel,iarr_wphys,iarr_sal,iarr_temp/)
tsrvars%ivarid = (/iarr_zeta,iarr_uvel,iarr_vvel,iarr_wphys,iarr_temp,iarr_airtemp/)
tsrvars%nrank = (/2,3,3,3,3,2/)

!2. Variable indices
!-------------------
!
ivarstsr(1:2,1) = (/1,6/)
ivarstsr(1:4,2) = (/2,3,4,5/)
!
!3. File parameters
!------------------
!

tsr2d(1)%defined = .TRUE.
!tsr3d(1)%defined = .FALSE.
!tsr2d(2:3)%defined = .FALSE.
tsr3d(2)%defined = .TRUE.

tsr2d(1)%filename = 'COHERENS_OUTPUT/best_bugHugo/best_2004_2004/'//TRIM(runtitle(:))//'_2D.nc'
tsr3d(2)%filename = 'COHERENS_OUTPUT/best_bugHugo/best_2004_2004/'//TRIM(runtitle(:))//'_3D.nc'

!4. Output grid
!--------------
tsrgpars(1)%time_format=1
tsrgpars(2)%time_format=1
tsrgpars(1)%nodim = 2
tsrgpars(2)%nodim = 3

tsrgpars(1)%tlims = (/0,nstep,60/)
!tsrgpars(1)%tlims = (/0,nstep,1440/)
tsrgpars(2)%tlims = (/0,nstep,60/) !1 hour output. se multiplica por timestep
!tsrgpars(2)%tlims = (/0,nstep,1440/) !1 hour output. se multiplica por timestep

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_tsr_params

!========================================================================

SUBROUTINE usrdef_tsr0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_tsr0d_vals* 0-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
!
! Description - test case bohai
!
! Reference -
!
! Calling program - time_series
!
! Module calls - define_out0d_vals
!
!************************************************************************
!
USE iopars
USE modids
USE model_output, ONLY: define_out0d_vals
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n0vars
REAL, INTENT(OUT), DIMENSION(n0vars) :: out0ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out0ddat* REAL    output data
!*n0vars*   INTEGER number of 0-D variables
!
!------------------------------------------------------------------------------
!

!PRINT *, 'second'
procname(pglev+1) = 'usrdef_tsr0d_vals'
CALL log_timer_in()


CALL define_out0d_vals(out0ddat,4,ivarid=(/iarr_ekin0d,iarr_epot0d,iarr_etot0d,&
                     & iarr_edissip0d/))

!---rescale to appropriate units
out0ddat(1:3) = 1.0E-09*out0ddat(1:3)
out0ddat(4) = 0.001*out0ddat(4)

CALL log_timer_out()   



RETURN

END SUBROUTINE usrdef_tsr0d_vals

!========================================================================

SUBROUTINE usrdef_tsr2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_tsr2d_vals* 2-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
!
! Description - test case bohai
!
! Reference -
!
! Calling program - time_series
!
! Module calls - define_out_2d_vals
!
!************************************************************************
!
USE modids
USE model_output, ONLY: define_out2d_vals

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, n2vars
REAL, INTENT(OUT), DIMENSION(n2vars) :: out2ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out2ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*n2vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!


CALL define_out2d_vals(out2ddat(3:4),i,j,2,ivarid=(/iarr_eflux2du,&
                     & iarr_eflux2dv/))
CALL define_out2d_vals(out2ddat(6:8),i,j,3,&
                     & ivarid=(/iarr_etot2d,iarr_edissip2d, iarr_sal/))
!PRINT *,'here2D'
!---rescale to appropriate units
out2ddat(3:4) = 1.0E-06*out2ddat(3:4)
out2ddat(6) = 1.0E-06*out2ddat(6)

RETURN

END SUBROUTINE usrdef_tsr2d_vals

!========================================================================

SUBROUTINE usrdef_tsr3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_tsr3d_vals* 3-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.1.2
!
! Description - test case bohai
!
! Reference -
!
! Calling program - time_series
!
! Module calls - define_out3d_vals
!
!************************************************************************
!
USE modids
USE sedids
USE switches
USE model_output, ONLY: define_out3d_vals
  
IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k, n3vars
REAL, INTENT(OUT), DIMENSION(n3vars) :: out3ddat



!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out3ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*k*        INTEGER Vertical index of output location
!*n3vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!

!---vector fields
!out3ddat(1) = Uvar_at_C(uvel(i:i+1,j,k),i,j,0,1)
!out3ddat(2) = Vvar_at_C(vvel(i,j:j+1,k),i,j,0,1)
!out3ddat(1) = define_out3d_vals(uvel(i:i+1,j,k),i,j,0,1)
!out3ddat(2) = define_out3d_vals(vvel(i,j:j+1,k),i,j,0,1)
!out3ddat(3) = wphys(i,j,k)

!---scalar field
!out3ddat(4) = sal(i,j,k)
!out3ddat(5) = ctot(i,j,k)
!CALL define_out3d_vals(out3ddat(8:10),i,j,k,8,ivarid=(/iarr_ctot,iarr_cvol,&
!                     & iarr_ctot,iarr_cvol/))
!---rescale to appropriate units
!out3ddat(4:7) = 1.0E-06*out3ddat(4:7)


RETURN

END SUBROUTINE usrdef_tsr3d_vals
