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

!=====================================================
module htm_v1_encoder
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v1_constants
 
    implicit none
    
        integer(INT32) :: time_step_max = 0
        logical(1), allocatable :: htm_encoder_data(:,:)
    
    contains
!=====================================================
    subroutine encode_pass_through(filename, dimension_1, dimension_2)
        character(len=*), intent(in) :: filename
        integer(INT32),   intent(in) :: dimension_1
        integer(INT32),   intent(in) :: dimension_2
        ! local variables
        integer(INT32) :: unitId, i1, i2, k, counter
        integer(INT8)  :: allocateStatus
        integer(INT8)  :: ios ! io success
        character(len=dimension_2+1) :: line
        open(newunit=unitId, file=filename, access='stream', form='formatted', action='read')

        !1] read the file without storing the content, but counting the time steps
        main_loop: do k = 1,100000
            do i1 = 1,dimension_1
                read(unit=unitId, fmt='(A)', iostat=ios), line
                if (ios == 0) then
                    !print '("INFO: encode_pass_through: time step ",i3,": line ",i4," has content ",A)', k, i1, line
                else 
                    time_step_max = k
                    !print '("INFO: encode_pass_through: done reading time step ",i5)', k
                    exit main_loop
                end if
            end do
            !print '("INFO: encode_pass_through: done reading time step ",i)', k
            ! read empy line
            !read(unit=unitId, fmt='(A)', iostat=ios), line
            !if (LEN(line) > 2) then
            !    print '("WARNING: encode_pass_through: expected an empty line but got line content ",A," with length=",i4)', line, LEN(line)
            !end if
        end do main_loop
        !print '("INFO: encode_pass_through: done reading ",i4," time steps")', time_step_max
        
        !2] create htm_encoder_data
        allocate(htm_encoder_data(dimension_1 * dimension_2, time_step_max), STAT = allocateStatus)
        if (allocateStatus /= 0) stop "*** encode_pass_through: Not enough memory ***"

        !3] reate the file again and store the content
        rewind(unit=unitId)
        do k = 1, time_step_max
            counter = 0
            do i1 = 1,dimension_1
                read(unit=unitId, fmt='(A)', iostat=ios), line
                if (ios == 0) then
                    !print '("INFO: encode_pass_through: read line ",i3," of time step ",i3," with counter=",A)', i1, k, line
                    do i2 = 1,dimension_2
                        counter = counter + 1
                        if (line(i2:i2) == '0') then
                            htm_encoder_data(counter, k) = .FALSE.
                        else if (line(i2:i2) == '1') then
                            htm_encoder_data(counter, k) = .TRUE.
                        else 
                            print '("WARNING: encode_pass_through: funky input at position ",i3," in line ",i4," with counter=",A)', i2, i1, line
                        end if
                    end do
                else
                    !ERROR
                end if
            end do
            ! read empy line
            read(unit=unitId, fmt='(A)', iostat=ios), line
        end do
        
        close(unit=unitId)
    end subroutine
!=====================================================
    subroutine fetch_sensor_data(t, sensor_data)
        integer(INT32), intent(in) :: t
        logical(1),     intent(out) :: sensor_data(:)
        integer(INT32) :: i
        i = MODULO(t-1, time_step_max-1)+1
        !print '("INFO: fetch_sensor_data: t=",i,"; time_step_max=",i,"; index=",i)',t, time_step_max, i
#       if _DEBUG
            if (.TRUE.) then
                print '("INFO:fetch_sensor_data; t=",i2,"; caseId=",i3)', t, i
                call print_input(i, SENSOR_DIM1, SENSOR_DIM2)
            end if
#       endif
        sensor_data(:) = htm_encoder_data(:,i)
    end subroutine
!=====================================================
    subroutine print_input(t, dimension_1, dimension_2)
        integer(INT32), intent(in) :: t
        integer(INT32),   intent(in) :: dimension_1
        integer(INT32),   intent(in) :: dimension_2

        integer(INT32) :: counter, i1, i2
    
        counter = 0
        do i1 = 1,dimension_1
            do i2 = 1,dimension_2
                counter = counter + 1
                if (htm_encoder_data(counter,t)) then
                    write(*, '("1")', advance='no')
                else 
                    write(*, '("0")', advance='no')
                end if
            end do
            write(*,*) ''
        end do
    end subroutine
!=====================================================
    subroutine encoder_cleanup()
        if (ALLOCATED(htm_encoder_data)) deallocate(htm_encoder_data)
    end subroutine encoder_cleanup
!=====================================================
end module htm_v1_encoder