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
module htm_v1_tp
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v1_constants
    use htm_v1_tools
 
    implicit none
    contains
!=====================================================
    subroutine temporal_pooling()
        integer(INT32) :: column_i
        do column_i = 1,N_COLUMNS
            call update_active_cell_state(column_i)
            call update_predictive_cell_state(column_i)
            call update_synapses(column_i)
        end do
    end subroutine temporal_pooling
!=====================================================
    subroutine update_active_cell_state(column)
        integer(INT32), intent(in) :: column
        ! local variables
        logical(1)     :: bottum_up_predicted, chosen
        integer(INT32) :: segment_i
        integer(INT8)  :: cell_i
        
        if (global_column_is_active(column)) then
            !print '("INFO: update_active_cell_state: column ",i4," is active." )', column
            bottum_up_predicted = .FALSE.
            chosen              = .FALSE.
            ! 1.1
            !do segment_i = 1, columns(column) % distal_dendrite_segment_size
            !    cell = columns(column) % distal_dendrite_segment_destination(segment_i)
            !    if (columns(column) % predictive_state_prev(cell)) then
            !    end if
            !end do
            
            do cell_i = 1,N_CELLS_PC
                if (columns(column) % predictive_state_prev(cell_i)) then
                    !print '("INFO: update_active_cell_state: column ",i4,"; cell ",i2," was predicted in the previous time step." )', column, cell_i

                    segment_i = get_active_sequence_segment_prev(cell_i, column) ! return (if present) an active sequence segment from cell.
                    !print '("INFO: update_active_cell_state: cell ",i2,"; column ",i4,"; get_active_sequence_segment_prev = ",i5 )', cell_i, column, segment_i

                    if (segment_i == -1) then
                        !print '("INFO: update_active_cell_state: column ",i4,"; cell ",i2," does not have an active sequence segment" )', column, cell_i
                    else 
                        !print '("INFO: update_active_cell_state: column ",i4,"; cell ",i2," was predicted and was a sequence segment" )', column, cell_i
                        bottum_up_predicted = .TRUE.
                        call set_cell_active(cell_i, column)

                        if (columns(column) % distal_dendrite_segment_active_prev(segment_i)) then
                            !print '("INFO: update_active_cell_state: column ",i4,"; cell ",i2," was predicted and was a sequence segment, segment was active" )', column, cell_i, segment_i
                            if (columns(column) % learn_state_prev(cell_i)) then
                                print '("INFO: update_active_cell_state: column ",i4,"; cell ",i2," was predicted and was a sequence segment, segment was active, cell is selected as learn_state_curr" )', column, cell_i, segment_i
                                chosen = .TRUE.
                                columns(column) % learn_state_curr(cell_i) = .TRUE. ! only one cell per column has learn state turned on, this cell matches the input best
                            end if
                        end if
                    end if
                end if
            end do
            ! 1.2
            if (.NOT.bottum_up_predicted) then
                !print '("INFO: update_active_cell_state: Bursting: no cells in column ",i4," have been predicted, making all active." )', column
                do cell_i = 1,N_CELLS_PC
                    call set_cell_active(cell_i, column) ! if no cells are in a predictive state, activate all cells in the column
                end do
            end if
            ! 1.3
            block
                integer(INT8)  :: winner_cell
                integer(INT32) :: best_segment, segment_match, segment_empty, synapse_i, cell2
                        
                if (.NOT.chosen) then
                    call get_best_matching_cell_and_segment_prev(column, winner_cell, best_segment, segment_match)
                    columns(column) % learn_state_curr(winner_cell) = .TRUE.

                    if (best_segment == -1) then
                        !print '("INFO: update_active_cell_state: column ",i4,": a new segment is added to cell ",i2,)', column, winner_cell
                        ! this column has no matching segment. Get_best_matching_cell_and_segment_prev returned the cell with the 
                        ! least number of segments. We add a new (random) segment to this cell.
                        call create_random_distral_dendrite_segment(winner_cell, column)
                    else
                        !print '("INFO: update_active_cell_state: column ",i4,"; segment ",i4," in cell ",i2," will be updated" )', column, best_segment, winner_cell
                        ! this columns has a matching segment, the best matching segment is part of cell winner_cell
                        
                        segment_empty = get_empty_distal_dendrite_segment_slot(column)
                        columns(column) % distal_dendrite_segment_tmp(segment_empty) % cell = winner_cell
                        columns(column) % distal_dendrite_segment_tmp(segment_empty) % is_sequence = .TRUE. ! the only place where is_sequence is written TRUE
                        columns(column) % distal_dendrite_segment_tmp(segment_empty) % idx = best_segment
                        columns(column) % distal_dendrite_segment_tmp(segment_empty) % origin(:)     = columns(column) % distal_dendrite_synapse_origin(:, best_segment)
                        columns(column) % distal_dendrite_segment_tmp(segment_empty) % permanance(:) = columns(column) % distal_dendrite_synapse_permanance(:, best_segment)

                        do synapse_i = 1, N_DD_SYNAPSES
                            cell2 = columns(column) % distal_dendrite_synapse_origin(synapse_i, best_segment)
                            if (.NOT.is_cell_active_prev(cell2)) then
                                columns(column) % distal_dendrite_segment_tmp(segment_empty) % origin(synapse_i)     = get_random_cell()
                                columns(column) % distal_dendrite_segment_tmp(segment_empty) % permanance(synapse_i) = INITIAL_PERMANENCE
                            end if
                        end do
                    end if
                end if
            end block
        end if
    end subroutine update_active_cell_state
!=====================================================
    subroutine get_best_matching_cell_and_segment_prev(column, best_cell, best_segment, segment_match)
        ! 
        integer(INT32), intent(in)  :: column
        integer(INT32), intent(out) :: best_segment, segment_match
        integer(INT8),  intent(out) :: best_cell
        ! local variables
        integer(INT32) :: segment_i, best_segment_tmp, n_active_synapses, n_active_synapses_tmp
        
        best_segment_tmp  = -1
        n_active_synapses_tmp = -1
        
        do segment_i = 1, columns(column) % distal_dendrite_segment_size
            n_active_synapses = columns(column) % distal_dendrite_segment_activity_prev(segment_i)
            if (n_active_synapses > n_active_synapses_tmp) then
                n_active_synapses_tmp = n_active_synapses
                best_segment_tmp      = segment_i
            end if
        end do
        
        !print '("INFO: get_best_matching_cell_and_segment_prev: the best matching segment for column ",i4," has match ",i4, "; threshold ",i5)', column, n_active_synapses_tmp, PROXIMAL_DENDRITE_SEGMENT_ACTIVATION_THRESHOLD
        
        if (n_active_synapses_tmp > PROXIMAL_DENDRITE_SEGMENT_ACTIVATION_THRESHOLD) then
            best_segment  = best_segment_tmp
            segment_match = n_active_synapses_tmp
            best_cell = columns(column) % distal_dendrite_segment_destination(best_segment_tmp)
        else
            ! if the best matching segment has a match lower then the segment activation threshold, then 
            ! return the cell with the fewest number of segments.
            best_segment  = -1
            segment_match = -1
            best_cell = get_cell_with_least_segments(column)
        end if
    end subroutine get_best_matching_cell_and_segment_prev
!=====================================================
    function get_cell_with_least_segments(column) result(r)
        integer(INT32), intent(in) :: column
        integer(INT8) :: r
        ! local variables
        integer(INT16) :: counter(N_CELLS_PC), selected_count
        integer(INT32) :: segment_i, cell_i
        integer(INT8)  :: selected_cell, cell
        
        counter = 0
        do segment_i = 1, columns(column) % distal_dendrite_segment_size
            cell = columns(column) % distal_dendrite_segment_destination(segment_i)
            counter(cell) = counter(cell) + 1
        end do
    
        selected_cell = 1
        selected_count = counter(1)
        do cell_i = 2, N_CELLS_PC
            if (counter(cell_i) < selected_count) then
                selected_count = counter(cell_i)
                selected_cell = cell_i
            end if
        end do
        r = selected_cell
    end function get_cell_with_least_segments
!=====================================================
    function get_active_sequence_segment_prev(cell, column) result(segment)
    ! return a segment index that is active and a sequence segment. If multiple segments are 
    ! applicable, return the sequence with the most activity are given preference. If no 
    ! sequences are applicable, return -1. 
        integer(INT8),  intent(in) :: cell
        integer(INT32), intent(in) :: column
        integer(INT32) :: segment
        ! local variable
        integer(INT16) :: activity, highest_activity
        integer(INT32) :: selected_segment, segment_i 


        !print '("get_active_sequence_segment_prev: entering: column ",i4,"; cell ",i2)', column, cell
        !call print_active_segments_prev(column)
        !call print_sequence_segments_prev(column)
        
        !somehow this method always return -1
        
        highest_activity = -1
        selected_segment = -1
        
        ! TODO: rather inefficient to loop through all segments and filter out only the segments for the provided cell
        do segment_i = 1,columns(column) % distal_dendrite_segment_size
            if (columns(column) % distal_dendrite_segment_destination(segment_i) == cell) then
                if (columns(column) % distal_dendrite_segment_active_prev(segment_i)) then
                    !print '("INFO: get_active_sequence_segment_prev: column ",i4," segment ",i4," was active in the previous time step.")', column, segment_i
                    if (columns(column) % distal_dendrite_segment_sequence_prev(segment_i)) then
                        
                        !print '("INFO: get_active_sequence_segment_prev: column ",i4,"; segment ",i4," was a sequence segment in the previous time step.")', column, segment_i
                       
                        activity = columns(column) % distal_dendrite_segment_activity_prev(segment_i)
                        if (activity > highest_activity) then
                            highest_activity = activity
                            selected_segment = segment_i
                        end if
                    end if
                end if
            end if
        end do
        segment = selected_segment
    end function get_active_sequence_segment_prev
!=====================================================
    subroutine update_predictive_cell_state(column)
        integer(INT32), intent(in) :: column
        integer(INT32) :: segment_i, synapse_i, segment_empty, best_segment
        integer(INT8)  :: cell
        integer(INT32) :: cell2
        real(REAL32)   :: old_permanance, new_permanance
        
        !print '("INFO: update_predictive_cell_state: entering")'
        call update_distal_dendrite_segment_activity(column)

        do segment_i = 1,columns(column) % distal_dendrite_segment_size
            !print '("INFO: update_predictive_cell_state: column ",i4," segment ",i6," has activity=",i4,"; threshold=",i4)', column, segment_i, columns(column) % distal_dendrite_segment_activity_curr(segment_i), DISTAL_DENDRITE_SEGMENT_ACTIVATION_THRESHOLD
            if (columns(column) % distal_dendrite_segment_active_curr(segment_i)) then
                !print '("INFO: update_predictive_cell_state: column ",i4," segment ",i6," is active (activity=",i4,"; threshold=",i4,")")', column, segment_i, columns(column) % distal_dendrite_segment_activity_curr(segment_i), DISTAL_DENDRITE_SEGMENT_ACTIVATION_THRESHOLD

                cell = columns(column) % distal_dendrite_segment_destination(segment_i)
                if (global_cell_is_active_curr(cell, column)) then 
                    !print '("INFO: update_predictive_cell_state: column ",i4," segment ",i6)', column, segment_i
                    columns(column) % predictive_state_curr(cell) = .TRUE. ! only place in which predictive state is written
                    
                    columns(column) % distal_dendrite_segment_sequence_curr(segment_i) = .TRUE.
                    
                    
                    if (LEARNING_ON) then
                        !----------------------------
                        !1] reinforce the currently active segment
                        segment_empty = get_empty_distal_dendrite_segment_slot(column)
                        columns(column) % distal_dendrite_segment_tmp(segment_empty) % cell = cell
                        columns(column) % distal_dendrite_segment_tmp(segment_empty) % is_sequence = .FALSE.
                        columns(column) % distal_dendrite_segment_tmp(segment_empty) % idx = segment_i
                        columns(column) % distal_dendrite_segment_tmp(segment_empty) % origin(:) = columns(column) % distal_dendrite_synapse_origin(:, segment_i)

                        do synapse_i = 1, N_DD_SYNAPSES
                            old_permanance = columns(column) % distal_dendrite_synapse_permanance(synapse_i, segment_i)
                            cell2 = columns(column) % distal_dendrite_synapse_origin(synapse_i, segment_i)
                            if (is_cell_active_curr(cell2)) then
                                new_permanance = MIN(1.0, old_permanance + PERMANENCE_INCREMENT)
                            else
                                new_permanance = MAX(0.0, old_permanance - PERMANENCE_DECREMENT)
                            end if
                            columns(column) % distal_dendrite_segment_tmp(segment_empty) % permanance(synapse_i) = new_permanance
                        end do
                        !----------------------------
                        !2] reinforce the segment that could have predicted this activation in the previous time step.
                        best_segment = get_best_matching_segment_prev(column, cell)
                        
                        if (best_segment == -1) then
                            !print '("INFO: update_predictive_cell_state: no best segment found.")'
                        else
                            !print '("INFO: update_predictive_cell_state: found a best segment ",i5,".")', best_segment
                            segment_empty = get_empty_distal_dendrite_segment_slot(column)
                            columns(column) % distal_dendrite_segment_tmp(segment_empty) % cell = cell
                            columns(column) % distal_dendrite_segment_tmp(segment_empty) % is_sequence = .FALSE.
                            columns(column) % distal_dendrite_segment_tmp(segment_empty) % idx = best_segment
                            columns(column) % distal_dendrite_segment_tmp(segment_empty) % origin(:)     = columns(column) % distal_dendrite_synapse_origin(:, best_segment)
                            columns(column) % distal_dendrite_segment_tmp(segment_empty) % permanance(:) = columns(column) % distal_dendrite_synapse_permanance(:, best_segment)

                            do synapse_i = 1, N_DD_SYNAPSES
                                cell2 = columns(column) % distal_dendrite_synapse_origin(synapse_i, segment_i)
                                if (.NOT.is_cell_active_prev(cell2)) then
                                    columns(column) % distal_dendrite_segment_tmp(segment_empty) % origin(synapse_i)     = get_random_cell()
                                    columns(column) % distal_dendrite_segment_tmp(segment_empty) % permanance(synapse_i) = INITIAL_PERMANENCE
                                end if
                            end do
                        end if
                    end if
                end if
            end if
        end do
    end subroutine update_predictive_cell_state
!=====================================================
    function get_best_matching_segment_prev(column, cell) result(best_segment)
        integer(INT32), intent(in) :: column
        integer(INT8),  intent(in) :: cell
        integer(INT32) :: best_segment
        ! local variables
        integer(INT32) :: segment_i, best_segment_tmp, n_active_synapses, n_active_synapses_tmp
        
        n_active_synapses_tmp = -1
        
        do segment_i = 1, columns(column) % distal_dendrite_segment_size
            if (columns(column) % distal_dendrite_segment_destination(segment_i) == cell) then
                n_active_synapses = columns(column) % distal_dendrite_segment_activity_prev(segment_i)
                if (n_active_synapses > n_active_synapses_tmp) then
                    n_active_synapses_tmp = n_active_synapses
                    best_segment_tmp      = segment_i
                end if
            end if
        end do
        
        if (n_active_synapses_tmp > PROXIMAL_DENDRITE_SEGMENT_ACTIVATION_THRESHOLD) then
            best_segment = best_segment_tmp
        else
            best_segment = -1
        end if
    end function get_best_matching_segment_prev
!=====================================================
    subroutine update_synapses(column)
        integer(INT32), intent(in) :: column
        integer(INT8)  :: cell
        integer(INT32) :: i
        
        do i = 1, N_DD_SEGMENTS_TMP_MAX
            cell = columns(column) % distal_dendrite_segment_tmp(i) % cell
            if (cell >= 0) then
                if (columns(column) % learn_state_curr(cell)) then
                    call adapt_distal_dendrite_segment(columns(column) % distal_dendrite_segment_tmp(i), .TRUE., column)
                    columns(column) % distal_dendrite_segment_tmp(i) % cell = -1
                else if (.NOT.columns(column) % predictive_state_curr(cell).AND.(columns(column) % predictive_state_prev(cell))) then
                    call adapt_distal_dendrite_segment(columns(column) % distal_dendrite_segment_tmp(i), .FALSE., column)
                    columns(column) % distal_dendrite_segment_tmp(i) % cell = -1
                else
                    !print '("INFO: update_synapses: leaving synapse ",i," untouched.")', i
                end if
            end if
        end do
    end subroutine update_synapses
!=====================================================
    subroutine adapt_distal_dendrite_segment(segment_tmp, positive_reinforcement, column)
        type(distal_dendrite_segment_tmp_t), intent(in) :: segment_tmp
        logical(1),     intent(in) :: positive_reinforcement
        integer(INT32), intent(in) :: column
        ! local variables
        integer(INT32) :: segment_i, synapse_i
        real(REAL32)   :: new_permanance
        integer(INT32) :: cell
        
        if (segment_tmp % idx < 0) then ! if idx is negative, segment_tmp is a new segment, not placed at a cell yet

            segment_i = columns(column) % distal_dendrite_segment_size
            if (segment_i < N_DD_SEGMENTS_MAX) then
                segment_i = segment_i + 1
                columns(column) % distal_dendrite_segment_size = segment_i
            else
                stop 'adapt_distal_dendrite_segment: too many segments'
            end if
            
            !print '("INFO: adapt_distal_dendrite_segment: column ",i4,", adding a new segment at position ",i5)', column, segment_i
            columns(column) % distal_dendrite_segment_destination(segment_i)  = segment_tmp % cell
            columns(column) % distal_dendrite_synapse_permanance(:,segment_i) = segment_tmp % permanance
            columns(column) % distal_dendrite_synapse_connected(:,segment_i) = (segment_tmp % permanance > CONNECTED_PERMANENCE)
            columns(column) % distal_dendrite_synapse_origin(:,segment_i)     = segment_tmp % origin
        else
            segment_i = segment_tmp % idx
            !print '("INFO: adapt_distal_dendrite_segment: column ",i4,", updating existing segment at pos ",i4,"; cell",i4)', column, segment_i, segment_tmp % cell
            
            columns(column) % distal_dendrite_segment_destination(segment_i)  = segment_tmp % cell
            columns(column) % distal_dendrite_synapse_permanance(:,segment_i) = segment_tmp % permanance
            columns(column) % distal_dendrite_synapse_connected(:,segment_i) = (segment_tmp % permanance > CONNECTED_PERMANENCE)
            columns(column) % distal_dendrite_synapse_origin(:,segment_i)     = segment_tmp % origin

            do synapse_i = 1,N_DD_SYNAPSES
                new_permanance = segment_tmp % permanance(synapse_i)
                cell = columns(column) % distal_dendrite_synapse_origin(synapse_i, segment_i)
                
                if (positive_reinforcement) then
                    !the currently synapses from currently active cells have contributed to correct prediction, and are awarded.
                    if (is_cell_active_curr(cell)) then
                        new_permanance = MIN(1.0, new_permanance + PERMANENCE_INCREMENT)
                    else
                        new_permanance = MAX(0.0, new_permanance - PERMANENCE_DECREMENT)
                    end if
                else 
                    if (is_cell_active_curr(cell)) then
                        new_permanance = MAX(0.0, new_permanance - PERMANENCE_DECREMENT)
                    else
                        new_permanance = MIN(1.0, new_permanance + PERMANENCE_INCREMENT)
                    end if
                end if
                columns(column) % distal_dendrite_synapse_permanance(synapse_i, segment_i) = new_permanance
            end do
            columns(column) % distal_dendrite_synapse_connected(:,segment_i) = (columns(column) % distal_dendrite_synapse_permanance(:,segment_i) > CONNECTED_PERMANENCE)
        end if
        columns(column) % distal_dendrite_segment_sequence_curr(segment_i) = segment_tmp % is_sequence
        
        if (.FALSE.) then
            if (columns(column) % distal_dendrite_segment_sequence_curr(segment_i)) then
                print '("INFO: adapt_distal_dendrite_segment: column ",i4,", segment ",i4," is a sequence segment")', column, segment_i
            end if
        end if
    end subroutine adapt_distal_dendrite_segment
!=====================================================
    subroutine create_random_distral_dendrite_segment(cell, column)
        integer(INT8),  intent(in) :: cell
        integer(INT32), intent(in) :: column
        integer(INT32) :: segment, synapse_i

        segment = get_empty_distal_dendrite_segment_slot(column)
        columns(column) % distal_dendrite_segment_tmp(segment) % cell = cell
        columns(column) % distal_dendrite_segment_tmp(segment) % is_sequence = .FALSE.
        columns(column) % distal_dendrite_segment_tmp(segment) % idx = -1

        do synapse_i = 1, N_DD_SYNAPSES
            columns(column) % distal_dendrite_segment_tmp(segment) % permanance(synapse_i) = INITIAL_PERMANENCE
            columns(column) % distal_dendrite_segment_tmp(segment) % origin(synapse_i)     = get_random_cell()
        end do
    end subroutine create_random_distral_dendrite_segment
!=====================================================
    function get_empty_distal_dendrite_segment_slot(column) result(r)
        integer(INT32), intent(in) :: column
        integer(INT32) :: r, segment_i
        r = -1
        do segment_i = 1, N_DD_SEGMENTS_TMP_MAX
            if (columns(column) % distal_dendrite_segment_tmp(segment_i) % cell < 0) then
                r = segment_i
                exit ! break from the do loop
            end if
        end do
        if (r == -1) then
            print '("WARNING: create_random_distral_dendrite_segment: could not find a empty slot for a new segment.")'
        end if
    end function
!=====================================================
    subroutine update_distal_dendrite_segment_activity(column)
        integer(INT32), intent(in) :: column
        ! local variables
        integer(INT32) :: segment_i, synapse_i, origin_cell
        integer(INT16) :: activity
        
        do segment_i = 1,columns(column) % distal_dendrite_segment_size
            activity = 0
            do synapse_i = 1,N_DD_SYNAPSES
                if (columns(column) % distal_dendrite_synapse_connected(synapse_i, segment_i)) then
                    !print '("INFO: update_distal_dendrite_segment_activity: column ",i4," segment ",i5,", synapse ",i5," is connected")', column, segment_i, synapse_i
                    
                    origin_cell = columns(column) % distal_dendrite_synapse_origin(synapse_i, segment_i)
                    if (is_cell_active_curr(origin_cell)) then
                        activity = activity + 1
                    end if
                end if
            end do
            
            !print '("INFO: update_distal_dendrite_segment_activity: column ",i4," segment ",i5,", activity=",i5)', column, segment_i, activity
            columns(column) % distal_dendrite_segment_activity_curr(segment_i) = activity
            columns(column) % distal_dendrite_segment_active_curr(segment_i) = (activity > DISTAL_DENDRITE_SEGMENT_ACTIVATION_THRESHOLD)
        end do
    end subroutine update_distal_dendrite_segment_activity
!=====================================================
end module htm_v1_tp
