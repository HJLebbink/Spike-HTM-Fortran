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
module htm_v2b_tools
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use IFPORT, only: IRAND, RAND
    use htm_v2b_constants
    implicit none
    contains
    
!========================================================================
    subroutine add_element_64(bitset64, pos)
        integer(INT64), intent(inout) :: bitset64
        integer(INT8),  intent(in)    :: pos
        bitset64 = IBSET(bitset64, pos)
    end subroutine add_element_64
!========================================================================
    subroutine clear_element_64(bitset64, pos)
        integer(INT64), intent(inout) :: bitset64
        integer(INT8),  intent(in)    :: pos
        bitset64 = IBCLR(bitset64, pos)
    end subroutine clear_element_64
!========================================================================
    elemental function get_element_64(bitset64, pos) result(r)
        integer(INT64), intent(in) :: bitset64
        integer(INT8),  intent(in) :: pos
        logical(1) :: r
        r = BTEST(bitset64, pos)
    end function get_element_64
!========================================================================
    function get_element_64_array(bitset, global_cell_id) result(r)
        integer(INT64), intent(in) :: bitset    (N_COLUMNS)
        type(global_cell_id_t), intent(in) :: global_cell_id
        logical(1) :: r
        r = BTEST(bitset(global_cell_id % column), global_cell_id % cell)
    end function get_element_64_array
!========================================================================
    subroutine clear_all_64(bitset64)
        integer(INT64), intent(out) :: bitset64
        bitset64 = 0
    end subroutine clear_all_64
!========================================================================
    subroutine clear_all_64_array(bitset64)
        integer(INT64), intent(out) :: bitset64(:)
        bitset64 = 0
    end subroutine clear_all_64_array
!========================================================================
    subroutine clear_all_logical(bitset)
        logical(1), intent(out) :: bitset(:)
        bitset = .FALSE.
    end subroutine clear_all_logical
!========================================================================
    subroutine set_all_64(bitset64)
        integer(INT64), intent(out) ::bitset64
        bitset64 = -1
    end subroutine set_all_64
!========================================================================
    subroutine add(bitset1, bitset2)
        logical(1), intent(in) :: bitset1(:)
        logical(1), intent(inout) :: bitset2(:)
        bitset2 = bitset2.OR.bitset1
    end subroutine add
!========================================================================
    subroutine add_element_logical(bitset, idx)
        logical(1),     intent(inout) :: bitset(:)
        integer(INT32), intent(in) :: idx
        bitset(idx) = .TRUE.
    end subroutine add_element_logical
!========================================================================
    subroutine clear_element_logical(bitset, idx)
        logical(1),     intent(inout) :: bitset(:)
        integer(INT32), intent(in) :: idx
        bitset(idx) = .FALSE.
    end subroutine clear_element_logical
!========================================================================
    function is_not_empty_64(bitset64) result(r)
        integer(INT64), intent(in) :: bitset64
        logical(1) :: r
        r = bitset64 /= 0
    end function is_not_empty_64
!========================================================================
    function is_not_empty_logical(bitset) result(r)
        logical(1), intent(in) :: bitset(:)
        logical(1) :: r
        r = NOT(ANY(bitset))
    end function is_not_empty_logical
!========================================================================
    function is_empty_logical(bitset) result(r)
        logical(1), intent(in) :: bitset(:)
        logical(1) :: r
        r = NOT(ANY(bitset))
    end function is_empty_logical
!========================================================================
    function is_empty_64(bitset64) result(r)
        integer(INT64), intent(in) :: bitset64
        logical(1) :: r
        r = bitset64 == 0
    end function is_empty_64
!========================================================================
    function is_empty_64_array(bitset64) result(r)
        integer(INT64), intent(in) :: bitset64(:)
        logical(1) :: r
        r = ALL(bitset64 == 0)
    end function is_empty_64_array
!========================================================================
    function count_active_64_array(bitset) result(r)
        integer(INT64), intent(in) :: bitset(:)
        integer(INT32) :: r, i
        r = 0
        do i = 1,SIZE(bitset)
            r = r + POPCNT(bitset(i))
        end do
    end function count_active_64_array
!========================================================================
    function count_active_logical(bitset) result(r)
        logical(1), intent(in) :: bitset(:)
        integer(INT32) :: r, i
        r = 0
        do i = 1,SIZE(bitset)
            if (bitset(i)) then
                r = r + 1
            end if
        end do
    end function count_active_logical
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
    function get_random_global_cell() result(r)
        type(global_cell_id_t) :: r
        r % column = rand_int32(1,N_COLUMNS)
        r % cell = rand_int32(1,N_CELLS_PC)
    end function get_random_global_cell
!========================================================================
    function get_random_global_cell_exclude(exclude_cells) result(global_cell_id)
        integer(INT64), intent(in) :: exclude_cells(N_COLUMNS)
        integer(INT32) :: try
        type(global_cell_id_t) :: global_cell_id
        
        do try = 1,10000
            global_cell_id = get_random_global_cell()
            if (NOT(get_element_64_array(exclude_cells, global_cell_id))) then
                exit
            end if
        end do
    end function get_random_global_cell_exclude
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
end module htm_v2b_tools
