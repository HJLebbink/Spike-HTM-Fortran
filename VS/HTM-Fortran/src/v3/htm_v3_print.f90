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
module htm_v3_print
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v3_constants
    use htm_v3_tools
    implicit none
    contains
    
!========================================================================
    subroutine print_bitset(bitset)
        logical(1), intent(in) :: bitset(:)
        integer(INT32) :: i
        write(*, '("(count ="(1X,I0)"):")', advance='no'), count_active(bitset)
        do i = 1,SIZE(bitset)
            if (bitset(i)) then
                write(*, '(1X,I0)', advance='no') i
!               write(*, '(i3)', advance='no') column_i
            end if
        end do
        write(*, '("")')
    end subroutine print_bitset
!========================================================================
    subroutine print_active_cells(bitset)
        logical(1), intent(in) :: bitset(N_CELLS_PC, N_COLUMNS)
        logical(1) :: active_column(N_COLUMNS)
        integer(INT32) :: column_i, n_active_columns, n_active_cells
        integer(INT8)  :: cell_i
        
        n_active_columns = 0
        do column_i = 1,N_COLUMNS
            active_column(column_i) = ANY(bitset(:,column_i))
            if (active_column(column_i)) then
                n_active_columns = n_active_columns + 1
            end if
        end do
        n_active_cells = count_active(RESHAPE(bitset, (/N_CELLS/)))
        
        write(*, '("INFO:number of active columns ="(1X,I0)", number of active cells ="(1X,I0)":")', advance='yes'), n_active_columns, n_active_cells
        if (n_active_columns > 0) then
            do column_i = 1,N_COLUMNS
                if (active_column(column_i)) then
                    write(*, '(i5)', advance='no'), column_i
                end if
            end do
            write(*, '("")')
            do cell_i = 1,N_CELLS_PC
                do column_i = 1,N_COLUMNS
                    if (active_column(column_i)) then
                        if (bitset(cell_i, column_i)) then
                            write(*, '("    X")', advance='no')
                        else
                            write(*, '("    .")', advance='no')
                        end if
                    end if
                end do
                write(*, '("")')
            end do
        end if
    end subroutine print_active_cells
!========================================================================
    subroutine print_active_columns(column_active, dimension_1, dimension_2)
        logical(1),     intent(in) :: column_active(N_COLUMNS)
        integer(INT32), intent(in) :: dimension_1
        integer(INT32), intent(in) :: dimension_2
        ! local variable
        integer(INT32) :: counter, i1, i2, sum
    
        sum = 0
        counter = 0
        do i1 = 1,dimension_1
            do i2 = 1,dimension_2
                counter = counter + 1
                if (counter <= N_COLUMNS) then
                    if (column_active(counter)) then
                        write(*, '("X")', advance='no')
                        sum = sum+1
                    else 
                        write(*, '(".")', advance='no')
                    end if
                else 
                    write(*, '("-")', advance='no')
                end if
            end do
            write(*,*) ''
        end do
        write(*, '("INFO number active columns"(1X,I0)", INHIBITION_TOP="(1X,I0))', advance='yes'), sum, INHIBITION_TOP
    end subroutine print_active_columns
!========================================================================
    subroutine print_sensor_activity(sensor_activity, dimension_1, dimension_2)
        logical(1),     intent(in) :: sensor_activity(N_SENSORS)
        integer(INT32), intent(in) :: dimension_1
        integer(INT32), intent(in) :: dimension_2
        ! local variable
        integer(INT32) :: counter, i1, i2
    
        counter = 0
        do i1 = 1,dimension_1
            do i2 = 1,dimension_2
                counter = counter + 1
                if (counter <= N_SENSORS) then
                    if (sensor_activity(counter)) then
                        write(*, '("X")', advance='no')
                    else 
                        write(*, '(".")', advance='no')
                    end if
                else
                    write(*, '("-")', advance='no')
                end if
            end do
            write(*,*) ''
        end do
    end subroutine print_sensor_activity
!========================================================================
    subroutine print_sensor_activity2(sensor_activity1, sensor_activity2, dimension_1, dimension_2)
        logical(1),     intent(in) :: sensor_activity1  (N_SENSORS)
        logical(1),     intent(in) :: sensor_activity2  (N_SENSORS)
        integer(INT32), intent(in) :: dimension_1
        integer(INT32), intent(in) :: dimension_2
        ! local variable
        integer(INT32) :: i1, i2
        
        logical(1) :: s1    (dimension_1, dimension_2)
        logical(1) :: s2    (dimension_1, dimension_2)

        s1 = RESHAPE(sensor_activity1, (/dimension_1, dimension_2/))
        s2 = RESHAPE(sensor_activity2, (/dimension_1, dimension_2/))
        
        do i2 = 1,dimension_2
            do i1 = 1,dimension_1
                if (s1(i1, i2)) then
                    write(*, '("X")', advance='no')
                else 
                    write(*, '(".")', advance='no')
                end if
            end do
            write(*, '(" | ")', advance='no')
            do i1 = 1,dimension_1
                if (s2(i1, i2)) then
                    write(*, '("X")', advance='no')
                else 
                    write(*, '(".")', advance='no')
                end if
            end do
            write(*, '(" | ")', advance='no')
            do i1 = 1,dimension_1
                if (s1(i1, i2) == s2(i1, i2)) then
                    write(*, '(".")', advance='no')
                else 
                    write(*, '("X")', advance='no')
                end if
            end do
            write(*,*) ''
        end do
    end subroutine print_sensor_activity2
!========================================================================
    subroutine print_inferred_sensor_activity(inferred_sensor_activity, dimension_1, dimension_2)
        integer(INT32), intent(in) :: inferred_sensor_activity(N_SENSORS)
        integer(INT32), intent(in) :: dimension_1
        integer(INT32), intent(in) :: dimension_2
        ! local variable
        integer(INT32) :: counter, i1, i2, sum
    
        counter = 0
        sum = 0
        do i1 = 1,dimension_1
            do i2 = 1,dimension_2
                counter = counter + 1
                write(*, '(i4)', advance='no') inferred_sensor_activity(counter)
                sum = sum + inferred_sensor_activity(counter)
            end do
            write(*,*) ''
        end do
        write(*,'("Sum="(1X,I0))'), sum
    end subroutine print_inferred_sensor_activity
!========================================================================
    subroutine print_projected_boost_factors(boost_factor, dimension_1, dimension_2)
        real(REAL32), intent(in) :: boost_factor    (N_SENSORS)
        integer(INT32), intent(in) :: dimension_1
        integer(INT32), intent(in) :: dimension_2
        ! local variable
        integer(INT32) :: counter, i1, i2
    
        counter = 0
        do i1 = 1,dimension_1
            do i2 = 1,dimension_2
                counter = counter + 1
                write(*, '(xf5.1)', advance='no') boost_factor(counter)
            end do
            write(*,*) ''
        end do
    end subroutine print_projected_boost_factors
!========================================================================
    subroutine print_boost_factors(boost_factor)
        real(REAL32), intent(in) :: boost_factor    (N_COLUMNS)
        ! local variable
        integer(INT32), parameter :: dimension_1 = 10
        integer(INT32) :: counter, i1

        counter = 0
        do while (counter < N_COLUMNS) 
            do i1 = 1,dimension_1
                counter = counter + 1
                if (counter <= N_COLUMNS) then
                    write(*, '(xf4.1)', advance='no') boost_factor(counter)
                end if
            end do
            write(*,*) ''
        end do
    end subroutine print_boost_factors
!========================================================================
    subroutine print_segment_activity(                                  &
            active_cells_all,                                           &
            permanence_threshold)
            
        logical(1),   intent(in) :: active_cells_all        (N_CELLS)
        real(REAL32), intent(in) :: permanence_threshold
        ! local variables
        integer(INT32) :: column_i, segment_i, synapse_i, n_active_synapses, global_cell_idx
        logical(1) :: something_to_print
        
        do column_i = 1, N_COLUMNS
            !1] check if somehting exists that we want to print
            something_to_print = .FALSE.
            do segment_i = 1, columns(column_i) % dd_segment_count
                do synapse_i = 1, N_DD_SYNAPSES
                    if (columns(column_i) % dd_synapse_permanence(synapse_i, segment_i) >= permanence_threshold) then
                        global_cell_idx = columns(column_i) % dd_synapse_origin(synapse_i, segment_i)
                        if (active_cells_all(global_cell_idx)) then
                            something_to_print = .TRUE.
                            exit
                        end if
                    end if
                end do
            end do
            
            if (something_to_print) then
                write(*, '("column"(1X,I0)": activity of"(1X,I0)" segments:")', advance='no'), column_i, columns(column_i) % dd_segment_count
                do segment_i = 1, columns(column_i) % dd_segment_count
                    n_active_synapses = 0
                    do synapse_i = 1, N_DD_SYNAPSES
                        if (columns(column_i) % dd_synapse_permanence(synapse_i, segment_i) >= permanence_threshold) then
                            global_cell_idx = columns(column_i) % dd_synapse_origin(synapse_i, segment_i)
                            if (active_cells_all(global_cell_idx)) then
                                n_active_synapses = n_active_synapses + 1
                            end if
                        end if
                    end do
                    write(*, '(1X,I0)', advance='no'), n_active_synapses
                end do
                write(*, '(".")', advance='yes')
            end if
        end do
    end subroutine print_segment_activity
end module htm_v3_print