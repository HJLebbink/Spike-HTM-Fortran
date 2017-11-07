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
module htm_v2_tools
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use IFPORT, only: IRAND, RAND
    use htm_v2_constants
    implicit none

        logical(1), parameter :: USE_LFSR_RANDOM_GENERATOR = .TRUE.
        integer(INT32)        :: current_random_number     = 'ACE1ACE1'X

    contains
!========================================================================
    ! get random integer equal (or larger) than min, and smaller (or equal) to max.
    function rand_int32(min, max) result(r)
        integer(INT32), intent(in) :: min, max
        integer(INT32) :: r

        if (USE_LFSR_RANDOM_GENERATOR) then
            r = rand_int32_lfsr(min, max)
        else
            r = rand_int32_internal(min, max)
        end if
#if _DEBUG
            if ((r < min).OR.(r > max)) print '("ERROR:tools:rand_int32: result"(1X,I0)"is invalid; min"(1X,I0)"; max"(1X,I0))', r, min, max
#endif
    
    end function rand_int32
!========================================================================
    function rand_int32_lfsr(min, max) result(r)
        integer(INT32), intent(in) :: min, max
        integer(INT32) :: r
        integer(INT64) :: u
        u = TRANSFER(current_random_number, INT64)
        r = min + MODULO(u, (max - min + 1))
        !if (max < 100) print '("INFO:tools:rand_int32: mod("(1X,I0)","(1X,I0)") ="(1X,I0))', u, (max - min + 1), r
        current_random_number = lfsr32_galois(current_random_number)
    end function rand_int32_lfsr
!========================================================================
    function rand_int32_internal(min, max) result(r)
        integer(INT32), intent(in) :: min, max
        integer(INT32) :: r
        real(REAL32) :: r2
        call RANDOM_NUMBER(r2)
        r = min + (r2 * (max - min + 1))
    end function rand_int32_internal
!========================================================================
    function rand_float(max) result(r)
        real(REAL32), intent(in) :: max
        real(REAL32) :: r
        if (USE_LFSR_RANDOM_GENERATOR) then
            block
                integer(INT64), parameter :: MAX_LFSR = INT(Z'7FFFFFFF', 8)
                r = (current_random_number * max) / MAX_LFSR;
#if _DEBUG
                    if (.FALSE.) print '("INFO:tools:rand_float: current_random_number ="(1X,I0)"; result ="(xf13.10))', current_random_number, r
#endif
                current_random_number = lfsr32_galois(current_random_number)
            end block
        else
            call RANDOM_NUMBER(r) ! Uniform distribution. The runtime-library implements the xorshift1024* random number generator (RNG). 
        end if
    end function rand_float
!========================================================================
    !Galois LFSR 32; see http://en.wikipedia.org/wiki/Linear_feedback_shift_register
    function lfsr32_galois(i) result(r)
        integer(INT32), intent(in) :: i
        integer(INT32) :: r
        ! local variables
        integer(INT32), parameter :: tabs = 'D0000001'X
        integer(INT32) :: j, lsb
        
        !return (i >> 1) ^ (-(signed int)(i & 1u) & 0xD0000001u); // The cast to signed int is needed to prevent the warning that unairy minus on unsigned ints yields unsigned ints.
        j = i
        lsb = -AND(j, 1)            ! Get LSB (i.e., the output bit).
        j = SHIFTR(j, 1)            ! Shift register
        j = XOR(j, AND(lsb, tabs));
        r = j
        
#if _DEBUG
        if (.FALSE.) print '("INFO:tools:lfsr32_galois: in"(1X,Z0)"; out"(1X,Z0))', i, j
#endif
    end function lfsr32_galois
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
!    function is_empty1(bitset) result(r)
!        logical(1), intent(in) :: bitset(:)
!        logical(1) :: r
!    this code has a bug
!        r = .NOT.(ANY(bitset))
!    end function is_empty1
!========================================================================
!    function is_empty2(bitset) result(r)
!        logical(1), intent(in) :: bitset(:,:)
!        logical(1) :: r
!    this code has a bug
!        r = .NOT.(ANY(bitset))
!    end function is_empty2
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
    function get_random_cell() result(r)
        integer(INT8) :: r
        r = INT(rand_int32(0,N_CELLS_PC-1)+1,INT8)
    end function get_random_cell
!========================================================================
    function get_random_global_cell() result(r)
        integer(INT32) :: r
        r = rand_int32(0,N_CELLS-1)+1
    end function get_random_global_cell
!========================================================================
    function get_random_global_cell_exclude(exclude_cells) result(r)
        logical(1), intent(in) :: exclude_cells(N_CELLS)
        integer(INT32) :: r
        integer(INT32) :: try
        r = -1
        do try = 1,10000
            r = get_random_global_cell()
            if (.NOT.(exclude_cells(r))) then
                return
            end if
        end do
        write(*, '("ERROR:get_random_global_cell_exclude: could not find random cell")')
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
        write(*, '("ERROR:get_random_global_cell_include: could not find random cell")')
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
end module htm_v2_tools
