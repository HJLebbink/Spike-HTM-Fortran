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
module htm_v3_tp2
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v3_constants
    use htm_v3_tools
    use htm_v3_print
    use htm_v3_tp_hist
    
    implicit none

        integer(INT32), parameter :: SYNAPSE_RECYCLE_T = 4000
    
    contains
!========================================================================
    subroutine compute_tp2(                                             &
            t,                                                          &
            active_columns,                                             &
            prev_active_cells_all,                                      &
            prev_winner_cells_all,                                      &
            prev_predictive_cells_all,                                  &
            enable_learn,                                               &
            ! out
            active_cells_all,                                           &
            winner_cells_all,                                           &
            predictive_cells_all)
            
        integer(INT32), intent(in) :: t
        logical(1), intent(in)  :: active_columns               (N_COLUMNS)
        logical(1), intent(in)  :: prev_active_cells_all        (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(in)  :: prev_winner_cells_all        (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(in)  :: prev_predictive_cells_all    (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(in)  :: enable_learn
        logical(1), intent(out) :: active_cells_all             (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: winner_cells_all             (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: predictive_cells_all         (N_CELLS_PC, N_COLUMNS)
        ! local variables
        logical(1)     :: predicted_columns              (N_COLUMNS)
        logical(1)     :: bursting_columns               (N_COLUMNS)
        integer(INT8)  :: cell
        integer(INT32) :: column_i, segment_i, segment
        integer(INT32) :: n_busting_columns, n_active_columns, n_active_synapses
        
        
#       if _DEBUG
            if (.FALSE.) then
                write(*, '("INFO:TP2:compute_tp2:prev_active_cells_all: ")', advance='no')
                call print_bitset(RESHAPE(prev_active_cells_all, (/N_CELLS/)))
            end if
            if (.FALSE.) then
                write(*, '("INFO:TP2:compute_tp2:prev_predictive_cells: ")', advance='no')
                call print_bitset(RESHAPE(prev_predictive_cells_all, (/N_CELLS/)))
            end if
#       endif
        
        !1] add the active_columns to the history
        !call prev_lrn_patterns_add(active_columns)
        
        !2] compute active cells
        n_busting_columns = 0
        n_active_columns  = 0
        bursting_columns(:) = .FALSE.
        
        do column_i = 1, N_COLUMNS
            predicted_columns(column_i) = ANY(prev_predictive_cells_all(:, column_i))
            if (active_columns(column_i)) then
                n_active_columns = n_active_columns + 1
                if (predicted_columns(column_i)) then
                    active_cells_all(:, column_i) = prev_predictive_cells_all(:, column_i)
                    winner_cells_all(:, column_i) = prev_predictive_cells_all(:, column_i)
                else
                    bursting_columns(column_i) = .TRUE.
                    n_busting_columns = n_busting_columns + 1
                    active_cells_all(:, column_i) = .TRUE. ! burst
                    
                    call best_matching_cell(column_i, prev_active_cells_all, cell)
                    winner_cells_all(:, column_i)    = .FALSE.
                    winner_cells_all(cell, column_i) = .TRUE.
                end if
            else
                active_cells_all(:, column_i) = .FALSE.
            end if
        end do

        if (2 * n_busting_columns > n_active_columns) then
            ! more than 50 per cent of active columns were not predicted!
            ! TODO do some backtracking
            
#           if _DEBUG
                if (.FALSE.) print '("INFO:TP2:compute_tp2: too many bursting columns. n_active_columns ="(1X,I0)"; n_bursting_columns ="(1X,I0))', n_active_columns, n_busting_columns
#           endif
        end if
        
#       if _DEBUG
            if (.TRUE.) then
                write(*, '("INFO:TP2:compute_tp2: active_cells_all =")', advance='yes')
                call print_active_cells(active_cells_all)
            end if
#       endif
        
        if (enable_learn) then
            ! update synapses
            
            ! if column bursts then
            !    Try to prevent the bursting from hapening agian. 
            !    1] Find a segment that almost predicted the previous active cells, 
            !       update this segment such that it may predict the previous cells.
            !       if no such segment can be found, add a new segment to a random 
            !       cell in this bursting column; this segment predicts the previous 
            !       active cells.
            ! else
            !    Update existing synapses to predict even better in the future.
            !    2] award synapses that predicted a cell in an active column,
            !    3] decay synapses that did not predict any cells,
            !    4] punish synapses that predicted a cell in an inactive column 
            ! end if
            
            do column_i = 1, N_COLUMNS
                
                columns(column_i) % prev_active_segments(1:columns(column_i) % dd_segment_count) = &
                    columns(column_i) % active_segments(1:columns(column_i) % dd_segment_count)
                
                if (bursting_columns(column_i)) then ! create new synapses
                    
                    call get_best_matching_segment(                                 &
                        column                  = column_i,                         &
                        active_cells_all        = prev_active_cells_all,            &
                        ! out
                        best_segment            = segment)
                    
                    if (segment == -1) then ! create new emtpy segment
                        call create_empty_DD_segment(                               &
                            column              = column_i,                         &
                            cell                = get_random_cell(),                &
                            ! out
                            segment             = segment)
                    end if
                    
                    if (segment /= -1) then
                        
                        call add_new_synapses(                                      &
                            column              = column_i,                         &
                            segment             = segment,                          &
                            winner_cells_all    = prev_winner_cells_all)
                        
                        columns(column_i) % dd_synapse_positive_activations(segment) = t
                        
                        if (.FALSE.) then
                            write(*, '("INFO:PT2:compute_tp2: column"(1X,I0)", prev_active_cells_all")', advance='yes'), column_i
                            call print_active_cells(prev_winner_cells_all)
                            !read *
                        end if
                    end if
                else ! update existing synapses
                    call update_synapses(                                           &
                        column                  = column_i,                         &
                        column_active           = active_columns(column_i),         &
                        column_predicted        = predicted_columns(column_i),      &
                        active_segments         = columns(column_i) % prev_active_segments,  &
                        prev_active_cells_all   = prev_active_cells_all)
                end if
            end do
        end if
        
        
        ! calculate the predicted cells 
        do column_i = 1, N_COLUMNS
            predictive_cells_all(:, column_i) = .FALSE.
            
            do segment_i = 1, columns(column_i) % dd_segment_count
                n_active_synapses = get_number_active_DD_synapses(          &
                    column                  = column_i,                     &
                    segment                 = segment_i,                    &
                    active_cells_all        = active_cells_all,             &
                    permanence_threshold    = DD_PERMANENCE_THRESHOLD)
#               if _DEBUG
                    if (.FALSE.) print '("INFO:TP2:compute_activity: column"(1X,I0)", segment"(1X,I0)": n_active_synapes ="(1X,I0)", DD_ACTIVATION_THRESHOLD ="(1X,I0)", DD_PERMANENCE_THRESHOLD =",f)', column_i, segment_i, n_active_synapses, DD_ACTIVATION_THRESHOLD, DD_PERMANENCE_THRESHOLD
#               endif
                if (n_active_synapses >= DD_ACTIVATION_THRESHOLD) then
                    cell = columns(column_i) % dd_segment_destination(segment_i)
                    predictive_cells_all(cell, column_i) = .TRUE.
                    columns(column_i) % active_segments(segment_i) = .TRUE.
                    
                    if (.FALSE.) then
                        if ((t - columns(column_i) % dd_synapse_last_active_iteration(segment_i)) > SYNAPSE_RECYCLE_T) then
                            print '("INFO:TP2:compute_activity: t"(1X,I0)": column"(1X,I0)", segment"(1X,I0)": had not been active for"(1X,I0)" time steps!")', t, column_i, segment_i, SYNAPSE_RECYCLE_T
                            !read *
                        end if
                    end if
                    
                    call inc(columns(column_i) % dd_synapse_total_activation(segment_i))
                    columns(column_i) % dd_synapse_last_active_iteration(segment_i) = t
                else
                    columns(column_i) % active_segments(segment_i) = .FALSE.
                endif
            end do
        end do
        
        ! analyse the use of segments
        if (.FALSE.) then
            block
                integer(INT32) :: synapse_used, synapse_i
                real(REAL32) :: total_num_segments, total_num_synapses
                total_num_segments = 0
                total_num_synapses = 0
                
                do column_i = 1, N_COLUMNS
                    
                    if (columns(column_i) % dd_segment_count > 0) then
                        !write(*, '("INFO:TP2:compute_tp2: t"(1X,I0)": column"(1X,I0)" has"(1X,I0)" segments")'), t, column_i, columns(column_i) % dd_segment_count
                        total_num_segments = total_num_segments + columns(column_i) % dd_segment_count
                    end if
                
                    do segment_i = 1, columns(column_i) % dd_segment_count
                        if ((t - columns(column_i) % dd_synapse_last_active_iteration(segment_i)) > SYNAPSE_RECYCLE_T) then
                            !write(*, '("INFO:TP2:compute_tp2: column"(1X,I0)", segment"(1X,I0)": at time"(1X,I0)" this segments has not been used for"(1X,I0)" time steps")'), column_i, segment_i, t, SYNAPSE_RECYCLE_T
                            call remove_segment(column_i, segment_i)
                        end if
                        
                        !write(*, '("INFO:TP2:compute_tp2: t"(1X,I0)": column"(1X,I0)", segment"(1X,I0)" has"(1X,I0)" synapses")'), t, column_i, segment_i, columns(column_i) % dd_synapse_count(segment_i)
                        total_num_synapses = total_num_synapses + columns(column_i) % dd_synapse_count(segment_i)
                        
                        synapse_used = 0
                        do synapse_i = 1, columns(column_i) % dd_synapse_count(segment_i)
                            if (columns(column_i) % dd_synapse_permanence(synapse_i, segment_i) > 0) then
                                call inc(synapse_used)
                            end if
                        end do
                        if (synapse_used < 5) then
                            !write(*, '("INFO:TP2:compute_tp2: t"(1X,I0)": column"(1X,I0)", segment"(1X,I0)" has only"(1X,I0)" used synapses")'), t, column_i, segment_i, synapse_used
                        end if
                    end do
                end do
                !write(*, '("INFO:TP2:compute_tp2: at time"(1X,I0)", average number of "f" segments, average number of "f" synapses per segment")'), t, (total_num_segments/N_COLUMNS), (total_num_synapses/total_num_segments)
            end block
        end if
        
#       if _DEBUG
            if (.TRUE.) then
                write(*, '("INFO:TP2:compute_tp2: predictive_cells_all =")', advance='yes')
                call print_active_cells(predictive_cells_all)
            end if
!            read *
#       endif
    end subroutine compute_tp2
!========================================================================
    function get_number_active_DD_synapses(                             &
            column,                                                     &
            segment,                                                    &
            active_cells_all,                                           &
            permanence_threshold) result(n_active_synapses)
    
        integer(INT32), intent(in) :: column, segment
        logical(1),     intent(in) :: active_cells_all          (N_CELLS)
        real(REAL32),   intent(in) :: permanence_threshold
        integer(INT32) :: n_active_synapses
        ! local variables
        integer(INT32) :: synapse_i, global_cell_id
        
        n_active_synapses = 0
        
#       if _DEBUG
        block
            integer(INT32) :: c1, c2
            c1 = columns(column) % dd_synapse_connected_count(segment)
            c2 = columns(column) % dd_synapse_count(segment)
            
            if ((c1 < 1).OR.(c1 > N_DD_SYNAPSES)) then
                print '("ERROR:TP2:get_number_active_DD_synapses: column"(1X,I0)"; segment"(1X,I0)"; dd_synapse_connected_count ="(1X,I0)" which is invallid")', column, segment, c1
                read *
            end if
            if ((c2 < 1).OR.(c2 > N_DD_SYNAPSES)) then
                print '("ERROR:TP2:get_number_active_DD_synapses: column"(1X,I0)"; segment"(1X,I0)"; dd_synapse_count ="(1X,I0)" which is invallid")', column, segment, c2
                read *
            end if
            if (c2 < c1) then
                print '("ERROR:TP2:get_number_active_DD_synapses: column"(1X,I0)"; segment"(1X,I0)"; dd_synapse_count ="(1X,I0)" and dd_synapse_connected_count="(1X,I0))', column, segment, c2, c1
                read *
            end if
            
            if (.TRUE.) print '("ERROR:TP2:get_number_active_DD_synapses: column"(1X,I0)"; segment"(1X,I0)"; dd_synapse_count ="(1X,I0)"; dd_synapse_connected_count="(1X,I0))', column, segment, c2, c1
            
            do synapse_i = 1, N_DD_SYNAPSES
                global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment)
                
                if (synapse_i <= c1) then  ! connected synapses
                    if (columns(column) % dd_synapse_permanence(synapse_i, segment) < permanence_threshold) then
                        print '("ERROR:TPR:get_number_active_DD_synapses: column"(1X,I0)"; segment"(1X,I0)"; synapse"(1X,I0)"; has permanence "f", which is below threshold "f". dd_synapse_connected_count="(1X,I0)".")', column, segment, synapse_i, columns(column) % dd_synapse_permanence(synapse_i, segment), permanence_threshold, c1
                        read *
                    end if
                    if ((global_cell_id < 1).OR.(global_cell_id > N_CELLS)) then
                        print '("ERROR:TP2:get_number_active_DD_synapses: column"(1X,I0)"; segment"(1X,I0)"; synapse"(1X,I0)"; has origin cell ="(1X,I0))', column, segment, synapse_i, global_cell_id
                        read *
                    end if
                else if (synapse_i <= c2) then  ! unconnected synapses
                    if (columns(column) % dd_synapse_permanence(synapse_i, segment) >= permanence_threshold) then
                        print '("ERROR:TPR:get_number_active_DD_synapses: column"(1X,I0)"; segment"(1X,I0)"; synapse"(1X,I0)"; has should have been connected but is was not")', column, segment, synapse_i
                        read *
                    end if
                    if ((global_cell_id < 1).OR.(global_cell_id > N_CELLS)) then
                        print '("ERROR:TP2:get_number_active_DD_synapses: column"(1X,I0)"; segment"(1X,I0)"; synapse"(1X,I0)"; has origin cell ="(1X,I0))', column, segment, synapse_i, global_cell_id
                        read *
                    end if
                else
                    ! unused synapses
                end if
            end do
        end block
#       endif
        
        do synapse_i = 1, columns(column) % dd_synapse_connected_count(segment)
            global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment)
            if (active_cells_all(global_cell_id)) then
                call inc(n_active_synapses)
            end if
        end do
        
#       if _DEBUG
            if (.FALSE.) print '("INFO:TP2:get_number_active_DD_synapses: column"(1X,I0)"; n_active_synapses"(1X,I0))', column, n_active_synapses
#       endif
    end function get_number_active_DD_synapses
!========================================================================
    subroutine update_synapses(                                         &
            column,                                                     &
            column_active,                                              &
            column_predicted,                                           &
            active_segments,                                            &
            prev_active_cells_all)
        
        integer(INT32), intent(in) :: column
        logical(1),     intent(in) :: column_active
        logical(1),     intent(in) :: column_predicted
        logical(1),     intent(in) :: active_segments       (N_DD_SEGMENTS_MAX)
        logical(1),     intent(in) :: prev_active_cells_all      (N_CELLS)
        ! local variables
        integer(INT32) :: segment_i, synapse_i, global_cell_id
        real(REAL32)   :: new_permanence, old_permanence
        
        if (column_predicted) then
            if (column_active) then
                do segment_i = 1, columns(column) % dd_segment_count
                    if (active_segments(segment_i)) then
                        do synapse_i = 1, columns(column) % dd_synapse_count(segment_i)
                            global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment_i)
                            old_permanence  = columns(column) % dd_synapse_permanence(synapse_i, segment_i)
                            if (prev_active_cells_all(global_cell_id)) then
                                new_permanence = inc_permanence(column, segment_i, synapse_i, old_permanence, DD_PERMANENCE_INC)
                            else
                                new_permanence = dec_permanence(column, segment_i, synapse_i, old_permanence, DD_PERMANENCE_DEC)
                            end if
                            columns(column) % dd_synapse_permanence(synapse_i, segment_i) = new_permanence
                        end do
                    end if
                end do
            else ! column is predicted but not active
                if (.TRUE.) then
                    do segment_i = 1, columns(column) % dd_segment_count
                        if (active_segments(segment_i)) then
                            do synapse_i = 1, columns(column) % dd_synapse_count(segment_i)
                                global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment_i)
                                if (prev_active_cells_all(global_cell_id)) then
                                    old_permanence = columns(column) % dd_synapse_permanence(synapse_i, segment_i)
                                    new_permanence = dec_permanence(column, segment_i, synapse_i, old_permanence, TP_DD_PREDICTED_SEGMENT_DEC)
                                    columns(column) % dd_synapse_permanence(synapse_i, segment_i) = new_permanence

#                                   if _DEBUG
                                        if (.FALSE.) print '("INFO:TP2:update_synapses: column"(1X,I0)"; segment"(1X,I0)"; synapse"(1X,I0)"; old_permanence =",f,"; new_permanence =",f)', column, segment_i, synapse_i, old_permanence, new_permanence
#                                   endif
                                end if
                            end do
                        end if
                    end do
                end if
            end if
        end if
    end subroutine
!========================================================================
    function inc_permanence(column, segment, synapse, permanence_in, increment) result(permanence_out)
        integer(INT32), intent(in) :: column, segment, synapse
        real(REAL32),   intent(in) :: permanence_in, increment
        real(REAL32) :: permanence_out
        
        permanence_out = MIN(1.0, permanence_in + increment)
        
        if ((permanence_in < DD_PERMANENCE_THRESHOLD).AND.(permanence_out >= DD_PERMANENCE_THRESHOLD)) then
            ! went from unconnected to connected!
            call inc(columns(column) % dd_synapse_connected_count(segment))
            call swap_synapses(column, segment, synapse, columns(column) % dd_synapse_connected_count(segment))
        end if
    end function
!========================================================================
    function dec_permanence(column, segment, synapse, permanence_in, decrement) result(permanence_out)
        integer(INT32), intent(in) :: column, segment, synapse
        real(REAL32),   intent(in) :: permanence_in, decrement
        real(REAL32) :: permanence_out
        integer(INT32) :: pos
        permanence_out = MAX(0.0, permanence_in - decrement)

        pos = columns(column) % dd_synapse_connected_count(segment)
        if (pos > 0) then
            if ((permanence_in >= DD_PERMANENCE_THRESHOLD).AND.(permanence_out < DD_PERMANENCE_THRESHOLD)) then
                ! went from connected to unconnected!
                call dec(columns(column) % dd_synapse_connected_count(segment))
                call swap_synapses(column, segment, synapse, pos)
            end if
        end if
    end function
!========================================================================
    subroutine swap_synapses(                                           &
            column, segment, synapse1, synapse2)
        integer(INT32), intent(in) :: column, segment, synapse1, synapse2
        
        integer(INT32) :: tmp1
        real(INT32)    :: tmp2
        
        tmp1 = columns(column) % dd_synapse_origin(synapse1, segment)
        columns(column) % dd_synapse_origin(synapse1, segment) = columns(column) % dd_synapse_origin(synapse2, segment)
        columns(column) % dd_synapse_origin(synapse2, segment) = tmp1

        tmp2 = columns(column) % dd_synapse_permanence(synapse1, segment)
        columns(column) % dd_synapse_permanence(synapse1, segment) = columns(column) % dd_synapse_permanence(synapse2, segment)
        columns(column) % dd_synapse_permanence(synapse2, segment) = tmp2
        
    end subroutine
!========================================================================
    subroutine add_new_synapses(                                        &
            column, segment,                                            &
            winner_cells_all)
    
        integer(INT32), intent(in)  :: column, segment
        logical(1),     intent(in)  :: winner_cells_all     (N_CELLS)
        ! local variables
        integer(INT32) :: i, synapse, n_synapses
        
        n_synapses = columns(column) % dd_synapse_count(segment)
        if ((n_synapses + NEW_NUMBER_SYNAPSES) > N_DD_SYNAPSES) then
            call cleanup_segment(column, segment)
        end if
        
        n_synapses = columns(column) % dd_synapse_count(segment)
        do i = 1, NEW_NUMBER_SYNAPSES
            if (n_synapses < N_DD_SYNAPSES) then
                call inc(n_synapses)
                columns(column) % dd_synapse_count(segment) = n_synapses
                
                synapse = n_synapses
                columns(column) % dd_synapse_permanence(synapse, segment) = DD_PERMANENCE_INIT
                columns(column) % dd_synapse_origin(synapse, segment) = &
                    get_random_global_cell_include(winner_cells_all)
                
                if (DD_PERMANENCE_INIT >= DD_PERMANENCE_THRESHOLD) then
                    call inc(columns(column) % dd_synapse_connected_count(segment))
                end if
            else
                if (.FALSE.) print '("INFO:TP2:add_new_synapses: column"(1X,I0)"; segment"(1X,I0)"; cannot add new synapses, segment if full")', column, segment
                return
            end if
        end do
    end subroutine
!========================================================================
    subroutine get_best_matching_segment(                               &
            column,                                                     &
            active_cells_all,                                           &
            best_segment)
    
        integer(INT32), intent(in)  :: column
        logical(1),     intent(in)  :: active_cells_all     (N_CELLS_PC, N_COLUMNS)
        integer(INT32), intent(out) :: best_segment
        ! local variables
        integer(INT32) :: n_segments, segment_i, n_active_synapses, best_activity
        
        best_segment = -1

        n_segments = columns(column) % dd_segment_count
        if (n_segments > 0) then
            best_activity = MIN_DD_ACTIVATION_THRESHOLD-1
            
            do segment_i = 1, n_segments
                n_active_synapses = get_number_active_DD_synapses(column, segment_i, active_cells_all, DD_PERMANENCE_THRESHOLD)
                if (n_active_synapses > best_activity) then
                    best_activity = n_active_synapses
                    best_segment  = segment_i
                end if
            end do
        end if
    end subroutine get_best_matching_segment
!========================================================================
    subroutine create_empty_DD_segment(                                 &
            column,                                                     &
            cell,                                                       &
            ! out
            segment)
    
        integer(INT32), intent(in)  :: column
        integer(INT8),  intent(in)  :: cell
        integer(INT32), intent(out) :: segment
        
        segment = columns(column) % dd_segment_count + 1
        if (segment > N_DD_SEGMENTS_MAX) then
            segment = -1
            print '("WARNING:TP2:create_empty_DD_segment: too many segments for column"(1X,I0)", ignoring new segment for cell"(1X,I0)".")', column, cell
        else
            columns(column) % dd_segment_destination(segment)   = cell
            columns(column) % dd_segment_count                  = segment

            columns(column) % dd_synapse_count(segment)         = 0
            columns(column) % dd_synapse_permanence(:, segment) = -1
            columns(column) % dd_synapse_origin(:, segment)     = -1
        end if
    end subroutine create_empty_DD_segment
!========================================================================
    subroutine best_matching_cell(                                      &
            column,                                                     &
            active_cells_all,                                           &
            ! out
            best_cell)

    ! get the best matching cell from the provided column given the active_cells
    ! if this code does not find a segment, then it is likely that a new segment 
    ! will be created
    
        integer(INT32), intent(in)  :: column
        logical(1),     intent(in)  :: active_cells_all     (N_CELLS)
        integer(INT8),  intent(out) :: best_cell
        ! local variables
        integer(INT32) :: segment_i, synapse_i
        integer(INT32) :: best_segment, num_active_synapses, highest_num_active_synapses, global_cell_id
  
        best_segment = -1
        highest_num_active_synapses = 0
        do segment_i = 1, columns(column) % dd_segment_count
            num_active_synapses = 0
            do synapse_i = 1, columns(column) % dd_synapse_count(segment_i)
                global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment_i)
                if (active_cells_all(global_cell_id)) then
                    num_active_synapses = num_active_synapses + 1
                end if
            end do
            if (num_active_synapses > highest_num_active_synapses) then
                highest_num_active_synapses = num_active_synapses
                best_segment  = segment_i
            end if
        end do

        if (best_segment /= -1) then
            best_cell = columns(column) % dd_segment_destination(best_segment)
        else
            best_cell = get_random_cell()
        end if
    end subroutine best_matching_cell
!========================================================================
    subroutine inc(i)
        integer(INT32), intent(inout) :: i
        i = i + 1
    end subroutine
!========================================================================
    subroutine dec(i) 
        integer(INT32), intent(inout) :: i
        i = i - 1
    end subroutine
!========================================================================
    subroutine remove_segment(                                          &
            column, segment)
        integer(INT32), intent(in) :: column, segment
        ! local variables
        integer(INT32) :: last
        
        last = columns(column) % dd_segment_count
        if (segment == last) then
            ! do nothing
            !print '("INFO:TP2:remove_segment: column"(1X,I0)", segment"(1X,I0)" is the last segment and is removed")', column, segment
        else
            !print '("INFO:TP2:remove_segment: column"(1X,I0)", segment"(1X,I0)" is removed and last"(1X,I0)" is copied")', column, segment, last
            columns(column) % dd_segment_destination(segment)   = columns(column) % dd_segment_destination(last)
            columns(column) % dd_synapse_permanence(:, segment) = columns(column) % dd_synapse_permanence(:, last)
            columns(column) % dd_synapse_origin(:, segment)     = columns(column) % dd_synapse_origin(:, last)
            columns(column) % dd_synapse_count(segment)         = columns(column) % dd_synapse_count(last)

            columns(column) % active_segments(segment)          = columns(column) % active_segments(last)
            columns(column) % prev_active_segments(segment)     = columns(column) % prev_active_segments(last)
            columns(column) % matching_segments(segment)        = columns(column) % matching_segments(last)
            columns(column) % prev_matching_segments(segment)   = columns(column) % prev_matching_segments(last)
            columns(column) % learning_segments(segment)        = columns(column) % learning_segments(last)
            columns(column) % sequence_segments(segment)        = columns(column) % sequence_segments(last)
            
            ! seg_updates are not used
        end if
        call dec(columns(column) % dd_segment_count)
        
        columns(column) % dd_synapse_total_activation(segment)      = columns(column) % dd_synapse_total_activation(last)
        columns(column) % dd_synapse_last_active_iteration(segment) = columns(column) % dd_synapse_last_active_iteration(last)
        columns(column) % dd_synapse_positive_activations(segment)  = columns(column) % dd_synapse_positive_activations(last)
    end subroutine
!========================================================================
    subroutine cleanup_segment(                                         &
            column, segment)
        integer(INT32), intent(in) :: column, segment
        ! local variables
        integer(INT32) :: synapse_i, i
        
        ! TODO: this is an inefficient code
        i = 0
        do synapse_i = 1, columns(column) % dd_synapse_count(segment)
            if (columns(column) % dd_synapse_permanence(synapse_i, segment) > 0) then
                i = i + 1
                columns(column) % dd_synapse_permanence(i, segment) = columns(column) % dd_synapse_permanence(synapse_i, segment)
                columns(column) % dd_synapse_origin(i, segment)     = columns(column) % dd_synapse_origin(synapse_i, segment)
                !print '("INFO:TP2:add_new_synapses: column"(1X,I0)"; segment"(1X,I0)"; synapse"(1X,I0)" has permanence == 0")', column, segment, synapse_i
            end if
        end do
        if (.FALSE.)  print '("INFO:TP2:cleanup_segment: column"(1X,I0)"; segment"(1X,I0)"; old number of segments"(1X,I0)", new number of segments"(1X,I0))', column, segment, columns(column) % dd_synapse_count(segment), i
        columns(column) % dd_synapse_count(segment) = i
        
    end subroutine cleanup_segment
!========================================================================
end module htm_v3_tp2
