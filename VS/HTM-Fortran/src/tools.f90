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
    module tools
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    implicit none
    contains
!=====================================================
    subroutine printElapsedTime(elapsedTimeCpu, elapsedTimeWall)
        real(REAL64),   intent(in) :: elapsedTimeCpu
        integer(INT64), intent(in) :: elapsedTimeWall
        integer(INT64) :: clockRate

        call SYSTEM_CLOCK(count_rate=clockRate)
        print '("CPU  time = ",f16.3," ms = ",f12.3," sec = ",f8.3," min = ",f6.3," hour")', &
                (elapsedTimeCpu)*1000, &
                (elapsedTimeCpu), &
                (elapsedTimeCpu)/60, &
                (elapsedTimeCpu)/(60*60)
        print '("Wall time = ",f16.3," ms = ",f12.3," sec = ",f8.3," min = ",f6.3," hour")', &
                (REAL(elapsedTimeWall,8)*1000)/clockRate, &
                (REAL(elapsedTimeWall,8))/clockRate, &
                (REAL(elapsedTimeWall,8)/60)/clockRate, &
                (REAL(elapsedTimeWall,8)/(60*60))/clockRate
    end subroutine
!=====================================================
    ! simple bubble sort
    pure subroutine sort(values)
        integer(INT32), intent(inout) :: values(:)
        ! local variables
        integer(INT32) tmp, i, j

        !DIR$ novector !Specifies that a particular loop should never be vectorized
        do i = 1, SIZE(values) - 1
            !DIR$ novector !Specifies that a particular loop should never be vectorized
            do j = i+1, SIZE(values)
                if (values(i).GT.values(j)) then
                    tmp = values(i)
                    values(i) = values(j)
                    values(j) = tmp
                end if
            end do
        end do
    end subroutine sort
!=====================================================
end module tools
