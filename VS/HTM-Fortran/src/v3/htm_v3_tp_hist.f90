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
module htm_v3_tp_hist
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v3_constants

        integer(INT32), parameter :: HISTORY_MAX_INF = 100
        integer(INT32), parameter :: HISTORY_MAX_LRN = 1000

        integer(INT32) :: history_size_inf = 0
        logical(1)     :: history_data_inf      (N_COLUMNS, HISTORY_MAX_INF)

        integer(INT32) :: history_size_lrn = 0
        logical(1)     :: history_data_lrn      (N_COLUMNS, HISTORY_MAX_LRN)

    ! NOTE: We don't use the same backtrack buffer for inference and learning
    ! because learning has a different metric for determining if an input from
    ! the past is potentially useful again for backtracking.
    !
    ! Our inference backtrack buffer. This keeps track of up to
    ! maxInfBacktrack of previous input. Each entry is a list of active column
    ! inputs.
    !self._prevInfPatterns = []

    ! Our learning backtrack buffer. This keeps track of up to maxLrnBacktrack
    ! of previous input. Each entry is a list of active column inputs
    !self._prevLrnPatterns = []

    contains
!========================================================================
    subroutine prev_inf_patterns_clear()
#       if _DEBUG
            if (.TRUE.) write(*, '("INFO:TP:prev_inf_patterns_clear")')
#       endif
        history_size_inf = 0
    end subroutine
!========================================================================
    subroutine prev_lrn_patterns_clear()
        history_size_lrn = 0
    end subroutine
!========================================================================
    function prev_inf_patterns_is_empty() result(r)
        logical(1) :: r
        r = (history_size_inf == 0)
    end function
!========================================================================
    function prev_lrn_patterns_is_empty() result(r)
        logical(1) :: r
        r = (history_size_lrn == 0)
    end function
!========================================================================
    function prev_inf_patterns_get_size() result(r)
        integer(INT32) :: r
        r = history_size_inf
    end function
!========================================================================
    function prev_lrn_patterns_get_size() result(r)
        integer(INT32) :: r
        r = history_size_lrn
    end function
!========================================================================
    function prev_inf_patterns_get_data(idx) result(r)
        integer(INT32), intent(in) :: idx
        logical(1) :: r(N_COLUMNS)
        r = history_data_inf(:, idx+1)
    end function
!========================================================================
    function prev_lrn_patterns_get_data(idx) result(r)
        integer(INT32), intent(in) :: idx
        logical(1) :: r(N_COLUMNS)
        r = history_data_lrn(:, idx+1)
    end function
!========================================================================
    subroutine prev_inf_patterns_remove_last()
#       if _DEBUG
            if (.TRUE.) write(*, '("INFO:TP:prev_inf_patterns_remove_last")')
#       endif
        if (history_size_inf > 0) then
            if (history_size_inf > 1) then ! only if more than 1 element exists we need to move data
                history_data_inf(:, 1:history_size_inf-1) = history_data_inf(:, 2:history_size_inf)
            end if
            history_size_inf = history_size_inf - 1
        end if
    end subroutine
!========================================================================
    subroutine prev_lrn_patterns_remove_last()
        if (history_size_lrn > 0) then
            if (history_size_lrn > 1) then ! only if more than 1 element exists we need to move data
                !TODO add a pointer to head, this to prevent copying
                history_data_lrn(:, 1:history_size_lrn-1) = history_data_lrn(:, 2:history_size_lrn)
            end if
            history_size_lrn = history_size_lrn - 1
        end if
    end subroutine
!========================================================================
    subroutine prev_inf_patterns_add(column_active) 
        logical(1), intent(in) :: column_active     (N_COLUMNS)
        history_size_inf = history_size_inf + 1
        if (history_size_inf <= HISTORY_MAX_INF) then
            history_data_inf(:, history_size_inf) = column_active
        else
            ! remove the last element
            call prev_inf_patterns_remove_last()
            history_data_inf(:, history_size_inf) = column_active
        end if
    end subroutine
!========================================================================
    subroutine prev_lrn_patterns_add(column_active) 
        logical(1), intent(in) :: column_active     (N_COLUMNS)
        history_size_lrn = history_size_lrn + 1
        if (history_size_lrn <= HISTORY_MAX_LRN) then
            history_data_lrn(:, history_size_lrn) = column_active
        else
            ! remove the last element
            call prev_lrn_patterns_remove_last()
            history_data_inf(:, history_size_inf) = column_active
        end if
    end subroutine
!========================================================================
end module htm_v3_tp_hist
