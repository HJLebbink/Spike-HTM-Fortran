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
    
module htm_v1_sp
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v1_constants
    use htm_v1_tools, only: print_active_columns
    implicit none
    contains
!=====================================================
    subroutine spatial_pooling()
        calc_overlap: block
            integer(INT32) :: column_i, synapse_i, overlap, sensor_idx
            
            do column_i = 1,N_COLUMNS
                overlap = 0
                do synapse_i = 1,N_PD_SYNAPSES
                    if (columns(column_i) % proximal_dendrite_synapse_permanance(synapse_i) > CONNECTED_PERMANENCE) then
                        sensor_idx = columns(column_i) % proximal_dendrite_synapse_origin(synapse_i)
                        if (global_sensor_activity(sensor_idx)) then
                            overlap = overlap + 1
                        end if
                    end if
                end do
                if (overlap < MIN_THRESHOLD) then
                    columns(column_i) % overlap = 0
                else 
                    columns(column_i) % overlap = overlap * columns(column_i) % boost
                end if
                !print '("INFO: spatial_pooling: overlap for column ",i4," is ",i)', column_i, columns(column_i) % overlap
            end do
        end block calc_overlap
        
        ! barrier: the overlap of all columns needs to be known
        
        inhibition_and_learning: block
            integer(INT32) :: column_i, cell_i, synapse_i, sensor_idx
            real(REAL32)   :: old_permanance, new_permanance
 
            do column_i = 1,N_COLUMNS
                global_column_is_active(column_i) = .FALSE.

                if (columns(column_i) % overlap > 0) then
                    if (column_has_most_overlap(column_i, MIN_THRESHOLD)) then
                        !print '("INFO: spatial_pooling: column ",i4," has most (overlap=",i4,").")', column_i, columns(column_i) % overlap
                        global_column_is_active(column_i) = .TRUE.

                        if (LEARNING_ON) then
                            do synapse_i = 1,N_PD_SYNAPSES
                                sensor_idx = columns(column_i) % proximal_dendrite_synapse_origin(synapse_i)
                                old_permanance = columns(column_i) % proximal_dendrite_synapse_permanance(synapse_i)

                                !print '("INFO: spatial_pooling:      old permanance for synpase ",i5," in column ",i4,":",f,").")', synapse_i, column_i, columns(column_i) % proximal_dendrite_synapse_permanance(synapse_i)
                                if (global_sensor_activity(sensor_idx)) then
                                    new_permanance = MIN(1.0, old_permanance + PERMANENCE_INCREMENT)
                                    !print '("INFO: spatial_pooling: inc: new permanance for synpase ",i5," in column ",i4,":",f,").")', synapse_i, column_i, permanance
                                else
                                    new_permanance = MAX(0.0, old_permanance - PERMANENCE_DECREMENT)
                                    !print '("INFO: spatial_pooling: dec: new permanance for synpase ",i5," in column ",i4,":",f,").")', synapse_i, column_i, permanance
                                end if
                                
                                if (.FALSE.) then
                                    if (old_permanance > CONNECTED_PERMANENCE) then
                                        if (new_permanance <= CONNECTED_PERMANENCE) then
                                            print '("INFO: spatial_pooling: column ",i3, ", synapse ",i4," went from connected to not connected")', column_i, synapse_i
                                        end if
                                    else
                                        if (new_permanance > CONNECTED_PERMANENCE) then
                                            print '("INFO: spatial_pooling: column ",i3, ", synapse ",i4," went from not connected to connected")', column_i, synapse_i
                                        end if
                                    end if
                                end if
                                columns(column_i) % proximal_dendrite_synapse_permanance(synapse_i) = new_permanance
                            end do
                        end if
                    end if
                end if
                
                ! TODO: 1] update network_boost(column_i)
                ! TODO: 2] increase permanence if input overlap is small
            end do
            
            ! TODO: recompute the inhibition radius
        end block inhibition_and_learning
        
        if (.FALSE.) then
            call print_active_columns()
            
        end if
        
    end subroutine spatial_pooling
!=====================================================
    function column_has_most_overlap(column, k) result(r)
        integer(INT32), intent(in) :: column, k
        logical(1) :: r
        ! local variables
        integer(INT32) :: n_columns_with_more_overlap, column_i, overlap_this_column
        
        overlap_this_column = columns(column) % overlap
        n_columns_with_more_overlap = 0
        
        if (GLOBAL_INHIBITION) then ! this can be done faster by just returing the k highest overlaps
            do column_i = 1,N_COLUMNS
                if (columns(column_i) % overlap > overlap_this_column) then
                    n_columns_with_more_overlap = n_columns_with_more_overlap + 1
                end if
            end do
        else
            !1] retrieve the columns in the neigborhood of column_i
            stop 'not implemented yet'
        end if
        !print '("INFO: column_has_most_overlap: column ",i4," has overlap=",i4,": ",i4," columns have more overlap (k=",i4,").")', column, overlap_this_column, n_columns_with_more_overlap, k
        r = (n_columns_with_more_overlap < k)
    end function column_has_most_overlap
!=====================================================
    subroutine update_inferred_sensor_activity()
        integer(INT32) :: column_i, synapse_i, input_idx
        global_inferred_sensor_activity = 0
        
        do column_i = 1,N_COLUMNS
            if (global_column_is_active(column_i)) then
                do synapse_i = 1, N_PD_SYNAPSES
                    if (columns(column_i) % proximal_dendrite_synapse_permanance(synapse_i) < PROXIMAL_DENDRITE_SEGMENT_ACTIVATION_THRESHOLD) then
                        input_idx = columns(column_i) % proximal_dendrite_synapse_origin(synapse_i)
                        global_inferred_sensor_activity(input_idx) = global_inferred_sensor_activity(input_idx) + 1
                    end if
                end do
            end if
        end do
        call print_inferred_sensor_activity()
    end subroutine
!=====================================================
    subroutine print_inferred_sensor_activity()
        integer(INT32) :: counter, i1, i2, sum
    
        counter = 0
        sum = 0
        do i1 = 1,SENSOR_DIM1
            do i2 = 1,SENSOR_DIM2
                counter = counter + 1
                write(*, '(i4)', advance='no') global_inferred_sensor_activity(counter)
                sum = sum + global_inferred_sensor_activity(counter)
            end do
            write(*,*) ''
        end do
        write(*,'("Sum=",i5)'), sum
    end subroutine
!=====================================================
end module htm_v1_sp
