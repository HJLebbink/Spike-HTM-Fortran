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
module htm_v3_tm
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v3_constants
    use htm_v3_tools
    use htm_v3_print
    implicit none

    contains
!========================================================================
    subroutine compute_tm(                                              &
            active_columns,                                             &
            prev_active_cells_all,                                      &
            prev_winner_cells_all,                                      &
            prev_predictive_cells_all,                                  &
            learn,                                                      &
            ! out
            active_cells_all,                                           &
            winner_cells_all,                                           &
            predictive_cells_all)
            
        logical(1), intent(in)  :: active_columns               (N_COLUMNS)
        logical(1), intent(in)  :: prev_active_cells_all        (N_CELLS)
        logical(1), intent(in)  :: prev_winner_cells_all        (N_CELLS)
        logical(1), intent(in)  :: prev_predictive_cells_all    (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(in)  :: learn
        logical(1), intent(out) :: active_cells_all             (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: winner_cells_all             (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: predictive_cells_all         (N_CELLS_PC, N_COLUMNS)
        
        ! local variables
        logical(1)     :: predicted_columns_local               (N_COLUMNS)
        integer(INT32) :: c_i, n_dd_segments
        
#       if _DEBUG
            if (.FALSE.) then
                write(*, '("INFO:TP:compute_tp:prev_predictive_cells: ")', advance='no')
                call print_bitset(RESHAPE(prev_predictive_cells_all, (/N_CELLS/)))
            end if
            if (.FALSE.) then
                write(*, '("INFO:TP:compute_tp:prev_winner_cells: ")', advance='no')
                call print_bitset(prev_winner_cells_all)
            end if
#       endif
        
        do c_i = 1, N_COLUMNS
            n_dd_segments = columns(c_i) % dd_segment_count
            columns(c_i) % prev_active_segments(1:n_dd_segments)    = columns(c_i) % active_segments(1:n_dd_segments)
            columns(c_i) % prev_matching_segments(1:n_dd_segments)  = columns(c_i) % matching_segments(1:n_dd_segments)
            columns(c_i) % prev_matching_cells     = columns(c_i) % matching_cells

            call activate_correctly_predictive_cells(                                   &
                prev_predictive_cells       = prev_predictive_cells_all(:, c_i),        &
                prev_matching_cells         = columns(c_i) % prev_matching_cells,       &
                active_column               = active_columns(c_i),                      &
                ! out
                active_cells                = active_cells_all(:, c_i),                 &
                winner_cells                = winner_cells_all(:, c_i),                 &
                predicted_column            = predicted_columns_local(c_i),             &
                predicted_inactive_cells    = columns(c_i) % predicted_inactive_cells)
        
            call burst_column(                                                          &
                column                      = c_i,                                      &
                active_column               = active_columns(c_i),                      &
                predicted_column            = predicted_columns_local(c_i),             &
                prev_active_cells_all       = prev_active_cells_all,                    &
                prev_winner_cells_all       = prev_winner_cells_all,                    &
                ! inout
                active_cells                = active_cells_all(:, c_i),                 &
                winner_cells                = winner_cells_all(:, c_i),                 &
                ! out
                learning_segments           = columns(c_i) % learning_segments)
            
            if (learn) then
                call learn_on_segments(                                                 &
                    column                   = c_i,                                     &
                    prev_active_segments     = columns(c_i) % prev_active_segments,     &
                    learning_segments        = columns(c_i) % learning_segments,        &
                    prev_active_cells_all    = prev_active_cells_all,                   &
                    winner_cells             = winner_cells_all(:, c_i),                &
                    prev_winner_cells_all    = prev_winner_cells_all,                   &
                    predicted_inactive_cells = columns(c_i) % predicted_inactive_cells, &
                    prev_matching_segments   = columns(c_i) % prev_matching_segments)
            end if
        end do
        
        ! barrier: active cells of all columns have to be known
        
        do c_i = 1, N_COLUMNS
            call compute_activity(                                                  &
                column                      = c_i,                                  &
                active_cells_all            = active_cells_all,                     &
                permanence_threshold        = DD_PERMANENCE_THRESHOLD,              &
                synapse_threshold           = DD_ACTIVATION_THRESHOLD,              &
                ! out
                active_segments             = columns(c_i) % active_segments,       &
                predictive_cells            = predictive_cells_all(:, c_i))

            call compute_activity(                                                  &
                column                      = c_i,                                  &
                active_cells_all            = active_cells_all,                     &
                permanence_threshold        = 0.0,                                  &
                synapse_threshold           = MIN_DD_ACTIVATION_THRESHOLD,          &
                ! out
                active_segments             = columns(c_i) % matching_segments,     &
                predictive_cells            = columns(c_i) % matching_cells)
        end do

#       if _DEBUG
            if (.FALSE.) then
                write(*, '("INFO:TP:compute_TP: number of predictive cells is"(1X,I0)"; synapse threshold is"(1X,I0)" synapse activity is:" )', advance='yes'), count_active(RESHAPE(predictive_cells_all, (/N_CELLS/))) , DD_ACTIVATION_THRESHOLD
                call print_segment_activity(active_cells_all, DD_PERMANENCE_THRESHOLD)
            end if
            if (.FALSE.) then
                print '("INFO:TP:compute_TP:predictive_cells:")'
                call print_bitset(RESHAPE(predictive_cells_all, (/N_CELLS/)))
            end if
#       endif
    end subroutine compute_tm
!========================================================================
    subroutine activate_correctly_predictive_cells(                     &
            prev_predictive_cells,                                      &
            prev_matching_cells,                                        &
            active_column,                                              &
            ! out
            active_cells,                                               &
            winner_cells,                                               &
            predicted_column,                                           &
            predicted_inactive_cells)

    ! Phase 1: activate the correctly predictive cells
    !   For the provided column, if the column is active then
    !       For all cells in the column, if the cell is predicted then
    !           mark the cell as active
    !           mark the cell as winner
    !           mark the column as predicted
    
        logical(1), intent(in)  :: prev_predictive_cells    (N_CELLS_PC)
        logical(1), intent(in)  :: prev_matching_cells      (N_CELLS_PC)
        logical(1), intent(in)  :: active_column
        logical(1), intent(out) :: active_cells             (N_CELLS_PC)
        logical(1), intent(out) :: winner_cells             (N_CELLS_PC)
        logical(1), intent(out) :: predicted_column
        logical(1), intent(out) :: predicted_inactive_cells (N_CELLS_PC)
    
#       if _DEBUG
            if (.FALSE.) then
                print '("INFO:TP:activate_correctly_predictive_cells: active column"(1X,I0)",; predicted column"(1X,I0)", prev predicted cells:")', predicted_column, is_not_empty1(prev_predictive_cells)
                call print_bitset(prev_predictive_cells)
            end if
#       endif
        
        if (active_column) then
            predicted_column = is_not_empty1(prev_predictive_cells)

#           if _DEBUG
                if (.FALSE.) print '("INFO:TP:activate_correctly_predictive_cells: active column; predicted"(1X,I0))', predicted_column
#           endif
            
            !do cell_i = 1,N_CELLS_PC
            !    if (prev_predictive_cells(cell_i)) then
            !        call add_element(active_cells, INT(cell_i, INT32))
            !        call add_element(winner_cells, INT(cell_i, INT32))
            !    end if
            !end do
            active_cells = prev_predictive_cells
            winner_cells = prev_predictive_cells
            call clear_all(predicted_inactive_cells)
        else
            predicted_column = .FALSE.

            if (TP_DD_PREDICTED_SEGMENT_DEC > 0.0) then
                ! predicted_inactive_cells are set such that these can be used
                ! in learn_on_segments to punish segments that make cells active
                ! while these cells were predicted to be inactive
                
                !do cell_i = 1,N_CELLS_PC
                !    if (prev_matching_cells(cell_i)) then
                !        call add_element(predicted_inactive_cells, INT(cell_i, INT32))
                !    end if
                !end do
                predicted_inactive_cells = prev_matching_cells
            end if
            
            call clear_all(active_cells)
            call clear_all(winner_cells)
        end if
    end subroutine activate_correctly_predictive_cells
!========================================================================
    subroutine burst_column(                                            &
            column,                                                     &
            active_column,                                              &
            predicted_column,                                           &
            prev_active_cells_all,                                      &
            prev_winner_cells_all,                                      &
            ! inout
            active_cells,                                               &
            winner_cells,                                               &
            ! out
            learning_segments)

        integer(INT32), intent(in)    :: column
        logical(1),     intent(in)    :: active_column
        logical(1),     intent(in)    :: predicted_column
        logical(1),     intent(in)    :: prev_active_cells_all  (N_CELLS)
        logical(1),     intent(in)    :: prev_winner_cells_all  (N_CELLS)
        logical(1),     intent(inout) :: active_cells           (N_CELLS_PC)
        logical(1),     intent(inout) :: winner_cells           (N_CELLS_PC)
        logical(1),     intent(out)   :: learning_segments      (N_DD_SEGMENTS_MAX)
        
!    public void burstColumns(ComputeCycle cycle, Connections c, 
!        Set<Column> activeColumns, 
!        Set<Column> predictedColumns, 
!        Set<Cell> prevActiveCells, 
!        Set<Cell> prevWinnerCells) 
!    {
!        activeColumns.removeAll(predictedColumns);
!        for(Column column : activeColumns) {
!            List<Cell> cells = column.getCells();
!            cycle.activeCells.addAll(cells);
!            CellSearch cellSearch = getBestMatchingCell(c, cells, prevActiveCells);
!            cycle.winnerCells.add(cellSearch.bestCell);
!            DistalDendrite bestSegment = cellSearch.bestSegment;
!            if(bestSegment == null && prevWinnerCells.size() > 0) {
!                bestSegment = cellSearch.bestCell.createSegment(c);
!            }
!            if(bestSegment != null) {
!                cycle.learningSegments.add(bestSegment);
!           }
!        }
!    }
        
        logical(1)     :: unpredicted_active_column, found_cell, found_segment
        integer(INT8)  :: best_cell
        integer(INT32) :: best_segment
 
        call clear_all(learning_segments(1:columns(column) % dd_segment_count))
        
        unpredicted_active_column = active_column.AND.NOT(predicted_column)
        if (unpredicted_active_column) then
#           if _DEBUG
                if (.FALSE.) print '("INFO:TP:burst_column: column",(1X,I0)," bursts")', column
#           endif
            
            call set_all(active_cells)  ! burst!
                
            call best_matching_cell(                                    &
                column                  = column,                       &
                active_cells_all        = prev_active_cells_all,        &
                ! out
                found_cell              = found_cell,                   &
                best_cell               = best_cell,                    &
                found_segment           = found_segment,                &
                best_segment            = best_segment)
                
            if (found_cell) then
                !print '("INFO:TP:burst_column: column",(1X,I0),"; cell",(1X,I0)," is added to winner cells")', column, best_cell
                call add_element(winner_cells, INT(best_cell, INT32)) ! only place that write to winner cells
            end if

            if (NOT(found_segment).AND.NOT(is_empty1(prev_winner_cells_all))) then
                
#               if _DEBUG
                    if (.FALSE.) then
                        write(*, '("INFO:TP:burst_column:prev_winner_cells: ")', advance='no')
                        call print_bitset(prev_winner_cells_all)
                    end if
#               endif
                call create_new_DD_segment(                             &
                    column                  = column,                   &
                    cell                    = best_cell,                &
                    prev_winner_cells_all   = prev_winner_cells_all)
                
            else
#               if _DEBUG
                    if (.FALSE.) print '("INFO:TP:burst_column: column",(1X,I0)," busts, but not adding a new segment")', column
#               endif
            end if
            
            if (found_segment) then
#               if _DEBUG
                    if (.FALSE.) print '("INFO:TP:burst_column: column",(1X,I0)," segment",(1X,I0)," is added to learning segments")', column, best_segment
#               endif
                call add_element(learning_segments, best_segment)
            end if
        end if
    end subroutine burst_column
!========================================================================
    subroutine learn_on_segments(                                       &
            column,                                                     &
            prev_active_segments,                                       &
            learning_segments,                                          &
            prev_active_cells_all,                                      &
            winner_cells,                                               &
            prev_winner_cells_all,                                      &
            predicted_inactive_cells,                                   &
            prev_matching_segments)

        integer(INT32), intent(in) :: column
        logical(1), intent(in) :: prev_active_segments      (N_DD_SEGMENTS_MAX)
        logical(1), intent(in) :: learning_segments         (N_DD_SEGMENTS_MAX)
        logical(1), intent(in) :: prev_active_cells_all     (N_CELLS)
        logical(1), intent(in) :: winner_cells              (N_CELLS_PC)
        logical(1), intent(in) :: prev_winner_cells_all     (N_CELLS)
        logical(1), intent(in) :: predicted_inactive_cells  (N_CELLS_PC)
        logical(1), intent(in) :: prev_matching_segments    (N_DD_SEGMENTS_MAX)
        
        logical(1)     :: is_learning_segment, is_from_winner_cell
        integer(INT32) :: segment_i
        integer(INT8)  :: cell
        
#      if _DEBUG
            if (.FALSE.) then
                if (is_not_empty1(prev_active_segments)) then
                    write(*, '("INFO:TP:learn_on_segments: column",(1X,I0),": prev_active_segments = ")', advance='no'), column
                    call print_bitset(prev_active_segments)
                end if
            end if
            if (.FALSE.) then
                if (is_not_empty1(learning_segments)) then
                    write(*, '("INFO:TP:learn_on_segments: column",(1X,I0),": learning_segments = ")', advance='no'), column
                    call print_bitset(learning_segments)
                end if
            end if
#      endif
        
        do segment_i = 1, columns(column) % dd_segment_count
            is_learning_segment = learning_segments(segment_i)
            if (prev_active_segments(segment_i).OR.(is_learning_segment)) then

                cell = columns(column) % dd_segment_destination(segment_i)
                is_from_winner_cell = winner_cells(cell)

#               if _DEBUG
                    if (.FALSE.) then
                        print '("INFO:TP:learn_on_segments: column",(1X,I0),"; segment",(1X,I0),": learning segment OR segment is previously active")', column, segment_i
                    end if
#               endif
                
                if (is_learning_segment.OR.is_from_winner_cell) then
                    call adapt_segment(                                     &
                        column                  = column,                   &
                        segment                 = segment_i,                &
                        prev_active_cells_all   = prev_active_cells_all,    &
                        permanence_inc          = DD_PERMANENCE_INC,        &
                        permanence_dec          = DD_PERMANENCE_DEC)
#                   if _DEBUG
                        if (.FALSE.) then
                            print '("INFO:TP:learn_on_segments: column",(1X,I0),"; cell",(1X,I0),"; segment",(1X,I0),": learning segment OR from winner cell")', column, cell, segment_i
                        end if
#                   endif
                end if
                if (is_learning_segment) then
                    call replace_DD_synapses(                               &
                        column                  = column,                   &
                        segment                 = segment_i,                &
                        permanence_threshold    = DD_PERMANENCE_THRESHOLD,  &
                        new_number_synapses     = NEW_NUMBER_SYNAPSES,      &
                        prev_winner_cells_all   = prev_winner_cells_all)
#                   if _DEBUG
                        if (.FALSE.) then
                            print '("INFO:TP:learn_on_segments: column",(1X,I0),"; cell",(1X,I0),"; segment",(1X,I0)," is a learning segment")', column, cell, segment_i
                        end if
#                   endif
                end if
            end if
        end do
        
        if (TP_DD_PREDICTED_SEGMENT_DEC > 0.0) then
            ! appreciate that this branch is taken for non active columns, thus we cannot parallelize on active columns
            
            ! punish the segments that predict the activity of a cell active 
            ! predicted_inactive_cells = prev_matching_cells
            
            do segment_i = 1, columns(column) % dd_segment_count
                if (prev_matching_segments(segment_i)) then
                    cell = columns(column) % dd_segment_destination(segment_i)
                    if (predicted_inactive_cells(cell)) then
                        call adapt_segment(                                         &
                            column                  = column,                       &
                            segment                 = segment_i,                    &
                            prev_active_cells_all   = prev_active_cells_all,        &
                            permanence_inc          = -TP_DD_PREDICTED_SEGMENT_DEC, &
                            permanence_dec          = 0.0)
#                       if _DEBUG
                            if (.FALSE.) then
                                if (TRACE_ON.AND.(column == TRACE_COLUMNS)) print '("INFO:TP:learn_on_segments: column",(1X,I0),"; segment",(1X,I0)," (from cell",(1X,I0),") is punished.")', column, segment_i, cell
                            end if
#                       endif
                    end if
                end if
            end do
        end if
    end subroutine learn_on_segments
!========================================================================
    subroutine compute_activity(                                        &
            column,                                                     &
            active_cells_all,                                           &
            permanence_threshold,                                       &
            synapse_threshold,                                          &
            ! out
            active_segments,                                            &
            predictive_cells)

    ! compute the active segments (and predictive cells) for the column, 
    ! given the provided active cells

        integer(INT32), intent(in)  :: column
        logical(1),     intent(in)  :: active_cells_all      (N_CELLS)
        real(REAL32),   intent(in)  :: permanence_threshold
        integer(INT32), intent(in)  :: synapse_threshold
        !out
        logical(1),     intent(out) :: active_segments       (N_DD_SEGMENTS_MAX)
        logical(1),     intent(out) :: predictive_cells      (N_CELLS_PC)

        ! local variables
        integer(INT8)  :: cell
        integer(INT32) :: segment_i, n_active_synapses
        
#       if _DEBUG
            if (.FALSE.) print '("INFO:TP:compute_activity: number of active cells =",(1X,I0))', count_active(active_cells_all)
#       endif
        
        call clear_all(predictive_cells)
        
        do segment_i = 1, columns(column) % dd_segment_count

            ! TODO: no need to count all active synapses, just count till synapse_threshold
            n_active_synapses = get_number_active_DD_synapses(          &
                column                  = column,                       &
                segment                 = segment_i,                    &
                active_cells_all        = active_cells_all,             &
                permanence_threshold    = permanence_threshold)
            
#           if _DEBUG
                if (.FALSE.) print '("INFO:TP:compute_activity: column",(1X,I0),", segment",(1X,I0),": n_active_synapes =",(1X,I0),", synapse_threshold =",(1X,I0),", permanence_threshold =",f)', column, segment_i, n_active_synapses, synapse_threshold, permanence_threshold
#           endif
            if (n_active_synapses >= synapse_threshold) then
                active_segments(segment_i) = .TRUE.
                cell = columns(column) % dd_segment_destination(segment_i)
                predictive_cells(cell) = .TRUE.
            else
                active_segments(segment_i) = .FALSE.
            end if
        end do
    end subroutine compute_activity
!========================================================================
    subroutine best_matching_cell(                                      &
            column,                                                     &
            active_cells_all,                                           &
            ! out
            found_cell,                                                 &
            best_cell,                                                  &
            found_segment,                                              &
            best_segment)

    ! get the best matching cell from the provided column given the active_cells
    ! if this code does not find a segment, then it is likely that a new segment 
    ! will be created
    
        integer(INT32), intent(in)  :: column
        logical(1),     intent(in)  :: active_cells_all     (N_CELLS)
        logical(1),     intent(out) :: found_cell
        integer(INT8),  intent(out) :: best_cell
        logical(1),     intent(out) :: found_segment
        integer(INT32), intent(out) :: best_segment
        ! local variables
        integer(INT32) :: segment_i, synapse_i, num_active_synapses, highest_num_active_synapses, global_cell_id
  
        found_cell    = .FALSE.
        found_segment = .FALSE.
        
        highest_num_active_synapses = 0
        do segment_i = 1, columns(column) % dd_segment_count
            num_active_synapses = 0
            do synapse_i = 1, N_DD_SYNAPSES
                global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment_i)
                if (global_cell_id == -1) exit
                if (active_cells_all(global_cell_id)) then
                    num_active_synapses = num_active_synapses + 1
                end if
            end do
            if (num_active_synapses > highest_num_active_synapses) then
                highest_num_active_synapses = num_active_synapses
                found_segment = .TRUE.
                best_segment  = segment_i
            end if
        end do

#       if _DEBUG
            if (.FALSE.) then
                if (found_segment) then
                    print '("INFO:TP:best_matching_cell: column",(1X,I0),"; found best segment",(1X,I0),".")', column, best_segment
                    !read *
                end if
            end if
#       endif
        
        if (found_segment) then
            best_cell     = columns(column) % dd_segment_destination(best_segment)
            found_cell    = .TRUE.
        end if
        
        if (NOT(found_cell)) then
            call least_used_cell(column, best_cell)
            found_cell    = .TRUE.
        end if
    end subroutine best_matching_cell
!========================================================================
    subroutine least_used_cell(                                         &
            column,                                                     &
            ! out
            best_cell)
    
        integer(INT32), intent(in)  :: column
        integer(INT8),  intent(out) :: best_cell
        ! local variables
        integer(INT16) :: counter               (N_CELLS_PC)
        integer(INT16) :: selected_count
        integer(INT32) :: segment_i
        integer(INT8)  :: cell_i
        
        !TODO: if multiple cells have an equal minimal usage, pick a random 
        ! cell from this set of cells 
        
        counter(:) = 0
        do segment_i = 1, columns(column) % dd_segment_count
            cell_i = columns(column) % dd_segment_destination(segment_i)
            counter(cell_i) = counter(cell_i) + 1
        end do
    
        best_cell = 1
        selected_count = counter(1)
        do cell_i = 2, N_CELLS_PC
            if (counter(cell_i) < selected_count) then
                selected_count = counter(cell_i)
                best_cell = cell_i
            end if
        end do
    end subroutine least_used_cell
!========================================================================
    subroutine adapt_segment(                                           &
            column,                                                     &
            segment,                                                    &
            prev_active_cells_all,                                      &
            permanence_inc,                                             &
            permanence_dec)
            
        integer(INT32), intent(in) :: column
        integer(INT32), intent(in) :: segment
        logical(1),     intent(in) :: prev_active_cells_all     (N_CELLS)
        real(REAL32),   intent(in) :: permanence_inc
        real(REAL32),   intent(in) :: permanence_dec
        
        integer(INT32) :: synapse_i, global_cell_id
        real(REAL32) :: new_permanence, old_permanence
        
        do synapse_i = 1, N_DD_SYNAPSES
            global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment)
            if (global_cell_id == -1) exit

            old_permanence  = columns(column) % dd_synapse_permanence(synapse_i, segment)
            if (prev_active_cells_all(global_cell_id)) then
                new_permanence = MIN(1.0, old_permanence + permanence_inc)
            else
                new_permanence = MAX(0.0, old_permanence - permanence_dec)
            end if
            columns(column) % dd_synapse_permanence(synapse_i, segment) = new_permanence
        end do
    end subroutine adapt_segment
!========================================================================
    subroutine replace_DD_synapses(                                     &
            column,                                                     &
            segment,                                                    &
            permanence_threshold,                                       &
            new_number_synapses,                                        &
            prev_winner_cells_all)
        
    ! original name: addNewSynapses. This code is not equal to the original 
    ! nupic code: up to NEW_NUMBER_SEGMENTS connections with a permanence 
    ! lower than the activation threshold are replaced with new (random) 
    ! connections.
        
        integer(INT32), intent(in) :: column
        integer(INT32), intent(in) :: segment
        real(REAL32),   intent(in) :: permanence_threshold
        integer(INT32), intent(in) :: new_number_synapses
        logical(1),     intent(in) :: prev_winner_cells_all     (N_CELLS)
        ! local variables
        integer(INT32) :: synapse_i, n_new_synapses, counter
        integer(INT32) :: selected_cells                     (N_DD_SYNAPSES)
        
        !1] count the number of synapses that will be replaced
        n_new_synapses = 0
        do synapse_i = 1, N_DD_SYNAPSES
            if (NOT(is_DD_connected(column, segment, synapse_i, permanence_threshold))) then
                n_new_synapses = n_new_synapses + 1
                if (n_new_synapses >= new_number_synapses) exit
            end if
        end do
        
        !2] get new cells from which we will get input
        call select_cell_to_learn_on(                           &
            column                  = column,                   &
            segment                 = segment,                  &
            prev_winner_cells_all   = prev_winner_cells_all,    &
            select_size             = n_new_synapses,           &
            ! out
            selected_cells       = selected_cells)

        !3] update the synapses with n_new_synapses new cells, 
        ! appreciate that n_new_synapses <= new_number_synapses
        counter = 0
        do synapse_i = 1, N_DD_SYNAPSES
            if (NOT(is_DD_connected(column, segment, synapse_i, permanence_threshold))) then
                counter = counter + 1
                columns(column) % dd_synapse_permanence(synapse_i, segment) = DD_PERMANENCE_INIT
                columns(column) % dd_synapse_origin(synapse_i, segment) = selected_cells(counter)

#               if _DEBUG
                    if (.FALSE.) print '("INFO:TP:replace_DD_synapses: column",(1X,I0),"; segment",(1X,I0),"; synapse",(1X,I0),"; counter",(1X,I0), "; new cell id",(1X,I0))', column, segment, synapse_i, counter, selected_cells(counter)
                    if ((selected_cells(counter) < 1).OR.(selected_cells(counter) > N_CELLS)) then
                        print '("ERROR:TP:replace_DD_synapses: column",(1X,I0),"; segment",(1X,I0),"; synapse",(1X,I0),"; counter",(1X,I0), "; new cell id",(1X,I0))', column, segment, synapse_i, counter, selected_cells(counter)
                        read *
                    end if
#               endif

                if (counter >= n_new_synapses) exit
            end if
        end do
        
#       if _DEBUG
            if (.FALSE.) print '("INFO:TP:replace_DD_synapses: column",(1X,I0),"; segment",(1X,I0),"; synapse",(1X,I0),"; done replacing",(1X,I0), " synapses")', column, segment, synapse_i, n_new_synapses
#       endif
            end subroutine replace_DD_synapses
!========================================================================
    subroutine create_new_DD_segment(                                   &
            column,                                                     &
            cell,                                                       &
            prev_winner_cells_all)
    
        integer(INT32), intent(in) :: column
        integer(INT8),  intent(in) :: cell
        logical(1),     intent(in) :: prev_winner_cells_all     (N_CELLS)
        ! local variable
        integer(INT32) :: new_segment_idx
        
        new_segment_idx = columns(column) % dd_segment_count + 1
        if (new_segment_idx > N_DD_SEGMENTS_MAX) then
            print '("WARNING:TP:create_new_DD_segment: too many segments for column"(1X,I0)", ignoring new segment for cell"(1X,I0)".")', column, cell
        else
#           if _DEBUG
                if (.FALSE.) print '("INFO:TP:create_new_DD_segment: column"(1X,I0)"; adding a new segment ("(1X,I0)") to cell"(1X,I0)"; total number of segments"(1X,I0))', column, new_segment_idx, cell, get_number_of_DD_segments(column, cell)
#           endif
            
            columns(column) % dd_segment_count = new_segment_idx
            
            columns(column) % dd_segment_destination(new_segment_idx) = cell
            columns(column) % dd_synapse_permanence(:, new_segment_idx) = -1
            columns(column) % dd_synapse_origin(:, new_segment_idx) = -1
        
            ! set all permanence to zero such that they are all replaced by add_new_synapses
            call replace_DD_synapses(                                       &
                column                  = column,                           &
                segment                 = new_segment_idx,                  &
                permanence_threshold    = DD_PERMANENCE_THRESHOLD,          &
                new_number_synapses     = 2*DD_ACTIVATION_THRESHOLD,        &
                prev_winner_cells_all   = prev_winner_cells_all)
                
        end if
    end subroutine create_new_DD_segment
!========================================================================
    subroutine select_cell_to_learn_on(                                 &
        column,                                                         &
        segment,                                                        &
        prev_winner_cells_all,                                          &
        select_size,                                                    &
        ! out
        selected_cells)

    ! pick global cell ids from winner_cells that do not already 
    ! have a pathway to segment in column, if not enough winner cells 
    ! exists, take non winner cells
    
        integer(INT32), intent(in)  :: column
        integer(INT32), intent(in)  :: segment
        logical(1),     intent(in)  :: prev_winner_cells_all    (N_CELLS)
        integer(INT32), intent(in)  :: select_size
        integer(INT32), intent(out) :: selected_cells        (N_DD_SYNAPSES)
        ! local variables

        integer(INT32) :: winner_cells_ptr(N_COLUMNS) ! assumption: should be large enough
        
        integer(INT32) :: winner_cells_ptr_i, synapse_i
        integer(INT32) :: rand_winner, rand_ptr, i, global_cell_id, counter
        logical(1)     :: already_present
        
        !1] make an array of winner cells
        winner_cells_ptr_i = 0
        do global_cell_id = 1, N_CELLS
            if (prev_winner_cells_all(global_cell_id)) then
                if (winner_cells_ptr_i < N_COLUMNS) then
                    winner_cells_ptr_i = winner_cells_ptr_i + 1
                    winner_cells_ptr(winner_cells_ptr_i) = global_cell_id
                else
                    print '("WARNING:TP:select_cell_to_learn_on: column",(1X,I0)," more winner cells than columns, ignoring thse winner cells")', column
                    exit
                end if
            end if
        end do

#       if _DEBUG
            if (.FALSE.) print '("INFO:TP:select_cell_to_learn_on: column",(1X,I0),"; found",(1X,I0)," previous winner cells")', column, winner_cells_ptr_i
#       endif
        
        counter = 0
        fill_selection: do i = 1, winner_cells_ptr_i
            rand_ptr = rand_int32(1, winner_cells_ptr_i)
            rand_winner                = winner_cells_ptr(rand_ptr)
            winner_cells_ptr(rand_ptr) = winner_cells_ptr(i) ! swap

            already_present = .FALSE.
            search_if_present: do synapse_i = 1, N_DD_SYNAPSES 
                if (columns(column) % dd_synapse_origin(synapse_i, segment) == rand_winner) then
                    already_present = .TRUE.
                    exit search_if_present
                end if
            end do search_if_present
            
            if (NOT(already_present)) then
                counter = counter + 1
                selected_cells(counter) = rand_winner
#               if _DEBUG
                    if (.FALSE.) print '("INFO:TP:select_cell_to_learn_on: column",(1X,I0),"; segment",(1X,I0),"; counter",(1X,I0), "; cell id",(1X,I0))', column, segment, counter, rand_winner
                    if ((rand_winner < 1).OR.(rand_winner > N_CELLS)) then
                        print '("ERROR:TP:select_cell_to_learn_on: column",(1X,I0),"; segment",(1X,I0),"; counter",(1X,I0), "; cell id",(1X,I0))', column, segment, counter, rand_winner
                    end if
#               endif
            end if

            if (counter >= select_size) exit fill_selection
        end do fill_selection
        
        !print '("INFO: TP: select_cell_to_learn_on: counter ",i4)', counter
        
        !3] not enough winner cells exists, add non-winner cells
        do i = counter+1, select_size
            selected_cells(i) = get_random_global_cell_exclude(prev_winner_cells_all) ! TODO duplicates can emerge
        end do
    end subroutine select_cell_to_learn_on
!========================================================================
    function get_number_active_DD_synapses(                             &
            column,                                                     &
            segment,                                                    &
            active_cells_all,                                           &
            permanence_threshold) result(r)
    
        integer(INT32), intent(in) :: column, segment
        logical(1),     intent(in) :: active_cells_all          (N_CELLS)
        real(REAL32),   intent(in) :: permanence_threshold
        integer(INT32) :: r
        ! local variables
        integer(INT32) :: synapse_i, global_cell_id, n_active_synapses
        
        n_active_synapses = 0
        
        do synapse_i = 1, N_DD_SYNAPSES
            global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment) ! 7% of time
            if (global_cell_id == -1) exit
            
            if (is_DD_connected(column, segment, synapse_i, permanence_threshold)) then ! 11% of time
#               if _DEBUG
                    if ((global_cell_id < 1).OR.(global_cell_id > N_CELLS)) then
                        print '("ERROR:TP:get_number_active_DD_synapses: column"(1X,I0)"; segment"(1X,I0)"; synapse"(1X,I0)"; has origin cell ="(1X,I0))', column, segment, synapse_i, global_cell_id
                    end if
#               endif
                if (active_cells_all(global_cell_id)) then ! this line accounts for 41% of time spend
                    n_active_synapses = n_active_synapses + 1
                end if
            end if
        end do

#       if _DEBUG
            if (.FALSE.) print '("INFO:TP:get_number_active_DD_synapses: column"(1X,I0)"; n_active_synapses"(1X,I0)"; N_DD_SYNAPSES"(1X,I0))', column, n_active_synapses, N_DD_SYNAPSES
#       endif
        r = n_active_synapses
    end function get_number_active_DD_synapses
!========================================================================
    function is_DD_connected(                                           &
            column,                                                     &
            segment,                                                    &
            synapse,                                                    &
            permanence_threshold) result(r)
    
        integer(INT32), intent(in) :: column, segment, synapse
        real(REAL32),   intent(in) :: permanence_threshold
        logical(1) :: r
        r = (columns(column) % dd_synapse_permanence(synapse, segment) >= permanence_threshold)
    end function is_DD_connected
!========================================================================
    function get_number_of_DD_segments(column, cell) result(r)
        ! INFO: no speed critical code, is only used in debug mode

        integer(INT32), intent(in) :: column
        integer(INT8),  intent(in) :: cell
        integer(INT32) :: r
        integer(INT32) :: segment_i
        r = 0
        do segment_i = 1, columns(column) % dd_segment_count
            if (cell == columns(column) % dd_segment_destination(segment_i)) then
                r = r + 1
            end if
        end do
    end function get_number_of_DD_segments
!========================================================================
    subroutine calc_anomaly_score(                                      &
            active_columns,                                             &
            prev_predicted_columns,                                     &
            anomaly_score)

    ! The raw anomaly score is the fraction of active columns not predicted.
    
        logical(1), intent(in) :: active_columns            (N_COLUMNS)
        logical(1), intent(in) :: prev_predicted_columns    (N_COLUMNS)
        real(REAL32), intent(out) :: anomaly_score          ! anomaly score 0..1 
        ! local variables
        integer(INT32) :: column_i, n_active_columns, n_active_and_predicted
        
        !TODO this could be done a lot faster with popcnt
        
        n_active_columns = 0
        n_active_and_predicted = 0
        
        do column_i = 1, N_COLUMNS
            if (active_columns(column_i)) then
                n_active_columns = n_active_columns + 1
                if (prev_predicted_columns(column_i)) then
                    n_active_and_predicted = n_active_and_predicted + 1
                end if
            end if
        end do
        
        if (n_active_columns > 0) then
            anomaly_score = REAL(n_active_columns - n_active_and_predicted) / REAL(n_active_columns)
        else
            if (ANY(prev_predicted_columns)) then
                anomaly_score = 1.0
            else 
                anomaly_score = 0.0
            end if
        end if
    end subroutine calc_anomaly_score
!========================================================================
end module htm_v3_tm