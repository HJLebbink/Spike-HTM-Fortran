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
module htm_v1_tools
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use IFPORT, only: IRAND, RAND
    use htm_v1_constants
 
    implicit none

    contains
!=====================================================
    function rand() result(r)
        real(REAL64) :: r
        call RANDOM_NUMBER(r) ! Uniform distribution
    end function rand
!=====================================================
    function rand_int32(min, max) result(r)
        integer(INT32) :: r
        integer(INT32), intent(in) :: min, max
        r = min + (rand() * ((max-min) + 1))
    end function rand_int32
!=====================================================
    function get_random_cell() result(r)
        integer(INT32) :: r
        integer(INT32) :: column, cell
        
        column = rand_int32(1,N_COLUMNS)
        cell   = rand_int32(1,N_CELLS_PC)
        r      = (column * N_CELLS_PC) + cell
    end function get_random_cell
!=====================================================
    function is_cell_active_curr(cell_id) result(r)
        integer(INT32), intent(in) :: cell_id
        logical(1) :: r
        integer(INT32) :: column, cell
        column = RSHIFT(cell_id, 6)+1 ! this only works N_CELLS_PC == 32
        cell = IAND(cell_id, 31)+1

        if (.FALSE.) then
            print '("is_cell_active_curr: cellId=",i6," => (column=",i4,"); active=",i1,".")', cell_id, column, cell
            call print_active_cells_curr()
        end if
        
        r = global_cell_is_active_curr(cell, column) 
        if (.FALSE..AND.r) then
            print '("is_cell_active_curr: cellId=",i6," => (column=",i4,"; cell=",i2,"); active=",i1,".")', cell_id, column, cell, r
        end if
    end function is_cell_active_curr
!=====================================================
    function is_cell_active_prev(cell_id) result(r)
        integer(INT32), intent(in) :: cell_id
        logical(1) :: r
        integer(INT32) :: column, cell
        column = RSHIFT(cell_id, 6)+1
        cell = IAND(cell_id, 31)+1
        !r = BTEST(active_cell_curr(column), cell) ! this only works N_CELLS_PC == 32
        r = global_cell_is_active_prev(cell, column) ! this only works N_CELLS_PC == 32
        if (.FALSE..AND.r) then
            print '("is_cell_active_prev: cellId=",i10," => (column=",i4,"; cell=",i2,"); active=",i1,".")', cell_id, column, cell, r
        end if
    end function is_cell_active_prev
!=====================================================
    subroutine set_cell_active_cell_id(cell_id) 
        integer(INT32), intent(in) :: cell_id
        integer(INT8)  :: cell
        integer(INT32) :: column
        column = RSHIFT(cell_id, 6)
        cell = IAND(cell_id, 31)
        call set_cell_active(cell, column)
    end subroutine set_cell_active_cell_id
!=====================================================
    subroutine set_cell_active(cell, column)
        integer(INT8),  intent(in) :: cell
        integer(INT32), intent(in) :: column
        !print '("set_cell_active: setting cell ",i2," in column ",i4," active")', cell, column
        global_cell_is_active_curr(cell, column) = .TRUE. ! only place where global_cell_is_active_curr is written
    end subroutine set_cell_active
!=====================================================
    subroutine to_string()
        integer(INT32) :: column_i
        do column_i = 1, N_COLUMNS
            call to_string_column(column_i)
        end do
    end subroutine to_string
!=====================================================
    subroutine to_string_column(column)
        integer(INT32), intent(in) :: column
        !integer(INT32) :: cell_i
    
        print '("INFO: to_string: column ",i4, " uses ",i4," distal dendrite segments." )', column, columns(column) % distal_dendrite_segment_size
        
        !do cell_i = 1, N_CELLS_PC
        !    print '("INFO: to_string: cell ",i2)', cell_i
        !    call to_string_cell(cell_i, column)
        !end do
    end subroutine to_string_column
!=====================================================
    subroutine to_string_cell(cell, column)
        integer(INT32), intent(in) :: cell, column
        integer(INT32) :: segment_i
    
        call to_string_proximal_dendrite_segment(column)
        do segment_i = 1, columns(column) % distal_dendrite_segment_size
            if (columns(column) % distal_dendrite_segment_destination(segment_i) == cell) then
                call to_string_distal_dendrite_segment(segment_i, column)
            end if
        end do
    end subroutine
!=====================================================
    subroutine to_string_proximal_dendrite_segment(column)
        integer(INT32), intent(in) :: column
    
        print '("INFO: PD synapse: origin", 10i)', columns(column) % proximal_dendrite_synapse_origin(:)
        print '("INFO: PD synapse: perman", 10f)', columns(column) % proximal_dendrite_synapse_permanance(:)
    end subroutine to_string_proximal_dendrite_segment
!=====================================================
    subroutine to_string_distal_dendrite_segment(segment, column)
        integer(INT32), intent(in) :: segment, column
    
        print '("INFO: DD synapse: origin", 10i)', columns(column) % distal_dendrite_synapse_origin(:,segment)
        print '("INFO: DD synapse: perman", 10f)', columns(column) % distal_dendrite_synapse_permanance(:,segment)
    end subroutine to_string_distal_dendrite_segment
!=====================================================
    subroutine print_active_columns()
        integer(INT32) :: column_i, n_active_columns
        n_active_columns = 0
        do column_i = 1,N_COLUMNS
            if (global_column_is_active(column_i)) then
                n_active_columns = n_active_columns + 1
            end if
        end do
        write(*, '("Active Columns (count="(1X,I0)"):")', advance='no'), n_active_columns
        do column_i = 1,N_COLUMNS
            if (global_column_is_active(column_i)) then
                write(*, '(1X,I0)', advance='no') column_i
!                 write(*, '(i3)', advance='no') column_i
            end if
        end do
        write(*, '(" ")')
    end subroutine
!=====================================================
    subroutine print_active_cells_curr()
        integer(INT32) :: column_i, cell_i
        
        do column_i = 1,N_COLUMNS
            do cell_i = 1,N_CELLS_PC
                if (global_cell_is_active_curr(cell_i, column_i)) then
                    write(*, '("c"(1X,I0)",cell"(1X,I0)")")', advance='no') column_i, cell_i
                end if
            end do
        end do
    end subroutine
!=====================================================
    subroutine print_sequence_segments_prev(column)
        integer(INT32), intent(in) :: column
        integer(INT32) :: segment_i
        integer(INT8)  :: cell
        write(*, '("print_sequence_segments_prev: column "(i4)":")', advance='yes'), column
        do segment_i = 1,N_DD_SYNAPSES
            if (columns(column) % distal_dendrite_segment_sequence_prev(segment_i)) then
                !write(*, '("(seg"(1X,I0)",cell"(1X,I0)")")', advance='no') segment_i, columns(column) % distal_dendrite_segment_destination(segment_i)
                cell = columns(column) % distal_dendrite_segment_destination(segment_i)
                print '("segment ",i4,"; cell ",i4)', segment_i, cell
            end if
        end do
        !write(*,*),''
    end subroutine
!=====================================================
    subroutine print_active_segments_prev(column)
        integer(INT32), intent(in) :: column
        integer(INT32) :: segment_i
        integer(INT8)  :: cell
        write(*, '("print_active_segments_prev: column "(i4)":")', advance='yes'), column
        do segment_i = 1,N_DD_SYNAPSES
            if (columns(column) % distal_dendrite_segment_active_prev(segment_i)) then
                !write(*, '("(seg"(1X,I0)",cell"(1X,I0)")")', advance='no') segment_i, columns(column) % distal_dendrite_segment_destination(segment_i)
                cell = columns(column) % distal_dendrite_segment_destination(segment_i)
                print '("segment ",i4,"; cell ",i4)', segment_i, cell
            end if
        end do
        !write(*,*),''
    end subroutine
!=====================================================
end module htm_v1_tools