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
module htm_v3_tools
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use IFPORT, only: IRAND, RAND
    use htm_v3_constants
    implicit none
    contains
!========================================================================
    subroutine clear_all(bitset)
        logical(1), intent(out) :: bitset(:)
        bitset = .FALSE.
    end subroutine clear_all
!========================================================================
    subroutine set_all(bitset)
        logical(1), intent(out) ::bitset(:)
        bitset = .TRUE.
    end subroutine set_all
!========================================================================
    subroutine add(bitset1, bitset2)
        logical(1), intent(in) :: bitset1(:)
        logical(1), intent(inout) :: bitset2(:)
        bitset2 = bitset2.OR.bitset1
    end subroutine add
!========================================================================
    subroutine add_element(bitset, idx)
        logical(1),     intent(inout) :: bitset(:)
        integer(INT32), intent(in) :: idx
        bitset(idx) = .TRUE.
    end subroutine add_element
!========================================================================
    subroutine clear_element(bitset, idx)
        logical(1),     intent(inout) :: bitset(:)
        integer(INT32), intent(in) :: idx
        bitset(idx) = .FALSE.
    end subroutine clear_element
!========================================================================
    function is_not_empty1(bitset) result(r)
        logical(1), intent(in) :: bitset(:)
        logical(1) :: r
        r = ANY(bitset)
    end function is_not_empty1
!========================================================================
    function is_empty1(bitset) result(r)
        logical(1), intent(in) :: bitset(:)
        logical(1) :: r
        r = NOT(ANY(bitset))
    end function is_empty1
!========================================================================
    function is_empty2(bitset) result(r)
        logical(1), intent(in) :: bitset(:,:)
        logical(1) :: r
        r = NOT(ANY(bitset))
    end function is_empty2
!========================================================================
    function count_active(bitset) result(r)
        logical(1), intent(in) :: bitset(:)
        integer(INT32) :: r, i
        r = 0
        do i = 1,SIZE(bitset)
            if (bitset(i)) then
                r = r + 1
            end if
        end do
    end function count_active
!========================================================================
    function rand() result(r)
        real(REAL64) :: r
        call RANDOM_NUMBER(r) ! Uniform distribution
    end function rand
!========================================================================
    function rand_int32(min, max) result(r)
        integer(INT32) :: r
        integer(INT32), intent(in) :: min, max
        r = min + (rand() * ((max-min) + 1))
    end function rand_int32
!========================================================================
    function get_random_cell() result(r)
        integer(INT8) :: r
        r = INT(rand_int32(1,N_CELLS_PC),INT8)
    end function get_random_cell
!========================================================================
    function get_random_global_cell() result(r)
        integer(INT32) :: r
        r = rand_int32(1,N_CELLS)
    end function get_random_global_cell
!========================================================================
    function get_random_global_cell_exclude(exclude_cells) result(r)
        logical(1), intent(in) :: exclude_cells(N_CELLS)
        integer(INT32) :: r
        integer(INT32) :: try
        r = -1
        do try = 1,10000
            r = get_random_global_cell()
            if (NOT(exclude_cells(r))) then
                return
            end if
        end do
        if (.TRUE.) write(*, '("ERROR:get_random_global_cell_exclude: could not find random cell")')
    end function get_random_global_cell_exclude
!========================================================================
    function get_random_global_cell_include(include_cells) result(r)
        logical(1), intent(in) :: include_cells(N_CELLS)
        integer(INT32) :: r
        integer(INT32) :: try
        r = -1
        do try = 1,10000
            r = get_random_global_cell()
            if (include_cells(r)) then
                return
            end if
        end do
        if (.FALSE.) write(*, '("ERROR:get_random_global_cell_include: could not find random cell")')
    end function get_random_global_cell_include
!========================================================================
    subroutine global_2_local_cell(global_cell_idx, cell, column)
        integer(INT32), intent(in)  :: global_cell_idx
        integer(INT8),  intent(out) :: cell
        integer(INT32), intent(out) :: column
        cell = IAND(global_cell_idx - 1, N_CELLS_PC-1) + 1 ! TODO: consider zero indexed arrays to get rid of the plus/minus 1
        column = RSHIFT(global_cell_idx - 1, N_BITS_CELL) + 1
    end subroutine global_2_local_cell
!========================================================================
    subroutine local_2_global_cell(cell, column, global_cell_idx)
        integer(INT8),  intent(in)  :: cell
        integer(INT32), intent(in)  :: column
        integer(INT32), intent(out) :: global_cell_idx
        global_cell_idx = IOR(LSHIFT(column - 1, N_BITS_CELL), IAND(cell - 1, N_CELLS_PC-1)) + 1
    end subroutine local_2_global_cell
!========================================================================
    subroutine test_global_cell()
        integer(INT32) :: column_i, global_cell_id
        integer(INT8)  :: cell_i
        
        do column_i = 1,10
            do cell_i = 1,N_CELLS_PC
                    call local_2_global_cell(cell_i, column_i, global_cell_id)
                    print '("INFO column ",i4,"; cell ",i2,"; global cell id ",i5)', column_i, cell_i, global_cell_id
            end do
        end do
    end subroutine 
end module htm_v3_tools
