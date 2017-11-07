! Fortran port of Nupic HTM
!
! Copyright (c) 2017 Henk-Jan Lebbink
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Affero Public License version 3 as
! published by the Free Software Foundation.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Affero Public License for more details.
!
! You should have received a copy of the GNU Affero Public License
! along with this program.  If not, see http://www.gnu.org/licenses.

!========================================================================
!  PROGRAM: htm_main
!  PURPOSE:  Entry point for the console application.
!========================================================================
program htm_main
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64

    use tools
    use htm_v2_constants
    use htm_v2_network

    implicit none

!========================================================================
    logical(1),         parameter :: progress       = .FALSE.
    integer(INT32),     parameter :: n_time_steps   = 400
    integer(INT32),     parameter :: n_times        = 50
!========================================================================
    real(REAL64)   :: timeStartCpu, timeEndCpu
    integer(INT64) :: timeStartWall, timeEndWall
    type(param_t)  :: param
!========================================================================
    ! code
    call CPU_TIME(timeStartCpu)
    call SYSTEM_CLOCK(timeStartWall)
    call RANDOM_SEED()
    
    call init_param(N_COLUMNS, N_BITS_CELL, param)
    print '("INFO: estimated memory usage"(1X,I0)" MB")', estimate_memory_usage_byte(param)/(1024*1024)
    
    call load_input()
    call allocate_memory(param)
    call run(n_time_steps, progress, n_times, param)
    
    call CPU_TIME(timeEndCpu)
    call SYSTEM_CLOCK(timeEndWall)
    
    print '("DONE: htm_main -----------------------------")'
    call printElapsedTime(timeEndCpu-timeStartCpu, timeEndWall-timeStartWall)
    print '("Press ENTER to close")'
    read *
    
end program htm_main
!========================================================================
