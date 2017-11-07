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
module htm_v3_tp
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v3_constants
    use htm_v3_tools
    use htm_v3_print
    use htm_v3_tp_hist
    
    implicit none

        integer(INT32), parameter :: MAX_SEQUENCE_LENGTH = 10
        
        real(REAL32) :: avg_input_density    = -1
    
        !internal_stats

        logical(1) :: g_inf_active_state          (N_CELLS_PC, N_COLUMNS)
        logical(1) :: g_inf_active_state_prev     (N_CELLS_PC, N_COLUMNS)
        
        logical(1) :: g_inf_predicted_state       (N_CELLS_PC, N_COLUMNS)
        logical(1) :: g_inf_predicted_state_prev  (N_CELLS_PC, N_COLUMNS)

        logical(1) :: g_lrn_active_state          (N_CELLS_PC, N_COLUMNS)
        logical(1) :: g_lrn_active_state_prev     (N_CELLS_PC, N_COLUMNS)
        
        logical(1) :: g_lrn_predicted_state       (N_CELLS_PC, N_COLUMNS)
        logical(1) :: g_lrn_predicted_state_prev  (N_CELLS_PC, N_COLUMNS)
        
        real(REAL32) :: g_cell_confidence         (N_CELLS_PC, N_COLUMNS)
        real(REAL32) :: g_cell_confidence_prev    (N_CELLS_PC, N_COLUMNS)
        
        real(REAL32) :: g_col_confidence          (N_COLUMNS)
        real(REAL32) :: g_col_confidence_prev     (N_COLUMNS)
        
        integer(INT32) :: learned_seq_length = 0
        integer(INT32) :: lrn_iteration_idx  = 0
        integer(INT32) :: iteration_idx      = 0
        real(REAL32)   :: avg_learned_seq_length = 0
        
        integer(INT32) :: pam_counter
        integer(INT32) :: pam_length
        
        logical(1) :: reset_called = .FALSE.
        
    contains
!========================================================================
    subroutine compute_tp(                                              &
            active_columns,                                             &
            enable_learn,                                               &
            collect_stats,                                              &
            compute_inf_output,                                         &
            ! out
            active_cells_all,                                           &
            predictive_cells_all)
            !tp_output)
            
        logical(1), intent(in)  :: active_columns               (N_COLUMNS)
        logical(1), intent(in)  :: enable_learn
        logical(1), intent(in)  :: collect_stats
        logical(1), intent(in)  :: compute_inf_output
        logical(1), intent(out) :: active_cells_all             (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: predictive_cells_all         (N_CELLS_PC, N_COLUMNS)
        !logical(1), intent(out) :: tp_output                    (N_CELLS_PC, N_COLUMNS)

        
        if (enable_learn) call inc(lrn_iteration_idx)
        call inc(iteration_idx)
        
        ! update the average input density
        if (avg_input_density == -1) then
            avg_input_density = count_active(active_columns)
        else
            avg_input_density = (0.99 * avg_input_density) + (0.01 * count_active(active_columns))
        end if
        
#       if _DEBUG
            if (.FALSE.) write(*, '("INFO:TP:compute_tp: avg_input_density",f)', advance='yes'),avg_input_density
#       endif
        
        if (compute_inf_output) then
            ! update past
            g_inf_active_state_prev           = g_inf_active_state
            g_inf_predicted_state_prev        = g_inf_predicted_state
            g_cell_confidence_prev            = g_cell_confidence
            g_col_confidence_prev             = g_col_confidence

            call update_inference_state(                                    &
                active_columns              = active_columns,               &
                cell_confidence_prev        = g_cell_confidence_prev,         &
                col_confidence_prev         = g_col_confidence_prev,          &
                inf_active_state_prev       = g_inf_active_state_prev,        &
                inf_predicted_state_prev    = g_inf_predicted_state_prev,     &
                ! out
                inf_active_state            = g_inf_active_state,             &
                inf_predicted_state         = g_inf_predicted_state,          &
                cell_confidence             = g_cell_confidence,              &
                col_confidence              = g_col_confidence)
        end if
        
        if (enable_learn) then
            ! update past
            g_lrn_predicted_state_prev = g_lrn_predicted_state
            g_lrn_active_state_prev    = g_lrn_active_state
            
            call update_learning_state(                                     &
                active_columns              = active_columns,               &
                lrn_predicted_state_prev    = g_lrn_predicted_state_prev,     &
                lrn_active_state_prev       = g_lrn_active_state_prev,        &
                ! out
                lrn_predicted_state         = g_lrn_predicted_state,          &
                lrn_active_state            = g_lrn_active_state)
        end if
    
        call apply_decay()
    
        ! update the prediction score stats
        if (collect_stats) then
            if (compute_inf_output) then
                !call update_stats_infer_end(                                &
                !    internal_stats      = internal_stats,                   &
                !    active_columns      = active_columns,                   &
                !    predicted_stats     = inf_predicted_state_prev,         &
                !    col_confidence      = col_confidence_prev)
            else
                !call update_stats_infer_end(                                &
                !    internal_stats      = internal_stats,                   &
                !    active_columns      = active_columns,                   &
                !    predicted_stats     = lrn_predicted_state_prev,         &
                !    col_confidence      = col_confidence_prev)
            end if
        end if
        
        
#       if _DEBUG
            if (.FALSE.) then
                write(*, '("INFO:TP:compute_tp: number of predictive cells is"(1X,I0)"; synapse threshold is"(1X,I0)" synapse activity is:" )', advance='yes'), count_active(RESHAPE(g_inf_predicted_state, (/N_CELLS/))) , DD_ACTIVATION_THRESHOLD
                call print_segment_activity(g_inf_predicted_state, DD_PERMANENCE_THRESHOLD)
            end if
            if (.FALSE.) then
                print '("INFO:TP:compute_tp:inf_active_state:")'
                call print_bitset(RESHAPE(g_inf_active_state, (/N_CELLS/)))
            end if
#       endif
        
        active_cells_all     = g_inf_active_state
        predictive_cells_all = g_inf_predicted_state
        !tp_output = inf_predicted_state.OR.inf_active_state
        reset_called = .FALSE.
    end subroutine compute_tp
!========================================================================
!========================================================================
    subroutine update_inference_state(                                  &
            active_columns,                                             &
            cell_confidence_prev,                                       &
            col_confidence_prev,                                        &
            inf_active_state_prev,                                      &
            
            inf_predicted_state_prev,                                   &
            ! out
            inf_active_state,                                           &
            inf_predicted_state,                                        &
            cell_confidence,                                            &
            col_confidence)

        logical(1), intent(in)    :: active_columns             (N_COLUMNS)
        logical(1), intent(in)    :: inf_active_state_prev      (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(inout)    :: inf_predicted_state_prev   (N_CELLS_PC, N_COLUMNS)
        real(REAL32), intent(in)  :: cell_confidence_prev       (N_CELLS_PC, N_COLUMNS)
        real(REAL32), intent(in)  :: col_confidence_prev        (N_COLUMNS)

        logical(1), intent(out)   :: inf_active_state           (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out)   :: inf_predicted_state        (N_CELLS_PC, N_COLUMNS)
        real(REAL32), intent(out) :: cell_confidence            (N_CELLS_PC, N_COLUMNS)
        real(REAL32), intent(out) :: col_confidence             (N_COLUMNS)
        
        ! local variables
        logical(1) :: in_sequence
        
        ! update our inference input history
        call prev_inf_patterns_add(active_columns)
        
        call infer_phase_1(                                             &
            active_columns              = active_columns,               &
            inf_predicted_state_prev    = inf_predicted_state_prev,     &
            use_start_cells             = reset_called,                 &
            ! out
            inf_active_state            = inf_active_state,             &
            in_sequence                 = in_sequence)
        
#       if _DEBUG
            if (.FALSE.) then
                write(*, '("INFO:TP:update_inference_state:inf_active_state: ")', advance='no')
                call print_active_cells(inf_active_state)
            end if
#       endif

        if (NOT(in_sequence)) then
            ! infer_backtrack will call infer_phase_2
            call infer_backtrack(                                       & 
                active_columns              = active_columns,           &
                ! inout
                inf_predicted_state_prev    = inf_predicted_state_prev, &
                inf_active_state            = inf_active_state,         &
                ! out
                inf_predicted_state         = inf_predicted_state,      &
                col_confidence              = col_confidence,           &
                cell_confidence             = cell_confidence)
        else 
            call infer_phase_2(                                         &
                inf_active_state            = inf_active_state,         &
                avg_input_density           = avg_input_density,        &
                ! out
                inf_predicted_state         = inf_predicted_state,      &
                col_confidence              = col_confidence,           &
                cell_confidence             = cell_confidence,          &
                in_sequence                 = in_sequence)
            
            if (NOT(in_sequence)) then
                call infer_backtrack(                                       & 
                    active_columns              = active_columns,           &
                    ! inout
                    inf_predicted_state_prev    = inf_predicted_state_prev, &
                    inf_active_state            = inf_active_state,         &
                    ! out
                    inf_predicted_state         = inf_predicted_state,      &
                    col_confidence              = col_confidence,           &
                    cell_confidence             = cell_confidence)
                end if
        end if
    end subroutine update_inference_state
!========================================================================
    subroutine infer_phase_1(                                           &
            active_columns,                                             &
            inf_predicted_state_prev,                                   &
            use_start_cells,                                            &
            ! out
            inf_active_state,                                           &
            in_sequence)
    
        logical(1), intent(in)  :: active_columns           (N_COLUMNS)
        logical(1), intent(in)  :: inf_predicted_state_prev (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(in)  :: use_start_cells
        logical(1), intent(out) :: inf_active_state         (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: in_sequence
        ! local variables
        integer(INT32) :: column_i, n_predicted_columns
        
        inf_active_state(:,:) = .FALSE.
        
        if (use_start_cells) then
            inf_active_state(1, :) = active_columns(:)
#           if _DEBUG
                if (.FALSE.) write(*, '("INFO:TP:infer_phase_1:entering, using start cells")')
#           endif
            in_sequence = .TRUE.
        else
            ! turn on any predicted cells in each column. If there are none,
            ! then turn on all cells (burst the column)
            n_predicted_columns = 0
            do column_i = 1, N_COLUMNS
                if (active_columns(column_i)) then
                    if (is_not_empty1(inf_predicted_state_prev(:, column_i))) then
                        inf_active_state(:, column_i) = inf_predicted_state_prev(:, column_i)
                        call inc(n_predicted_columns)
                    else
                        inf_active_state(:, column_i) = .TRUE.  ! whole column bursts
                    end if
                end if
            end do
#           if _DEBUG
                if (.TRUE.) write(*, '("INFO:TP:infer_phase_1:entering, not using start cells: n_predicted_columns"(1X,I0))'), n_predicted_columns
                !read *
#           endif
            ! did we predict this input well enough?
            in_sequence = (n_predicted_columns >= (0.5 * count_active(active_columns)))
        end if
    end subroutine infer_phase_1
!========================================================================
    subroutine infer_phase_2(                                           &
            inf_active_state,                                           &
            avg_input_density,                                          &
            ! out
            inf_predicted_state,                                        &
            col_confidence,                                             &
            cell_confidence,                                            &
            in_sequence)
    
    !Phase 2 for the inference state. The computes the predicted state, then
    !checks to insure that the predicted state is not over-saturated, i.e.
    !look too close like a burst. This indicates that there were so many
    !separate paths learned from the current input columns to the predicted
    !input columns that bursting on the current input columns is most likely
    !generated mix and match errors on cells in the predicted columns. If
    !we detect this situation, we instead turn on only the start cells in the
    !current active columns and re-generate the predicted state from those.

    !@returns True if we have a decent guess as to the next input.
    !         Returing False from here indicates to the caller that we have
    !         reached the end of a learned sequence.

        logical(1),   intent(in)  :: inf_active_state       (N_CELLS_PC, N_COLUMNS)
        real(REAL32), intent(in)  :: avg_input_density
        
        logical(1),   intent(out) :: inf_predicted_state    (N_CELLS_PC, N_COLUMNS)
        real(REAL32), intent(out) :: col_confidence         (N_COLUMNS)
        real(REAL32), intent(out) :: cell_confidence        (N_CELLS_PC, N_COLUMNS)
        logical(1),   intent(out) :: in_sequence
        
        ! local variables
        integer(INT32) :: column_i, segment_i, dc
        integer(INT8)  :: cell
        integer(INT32) :: n_predicted_columns, n_active_synapses
        real(REAL32)   :: sum_confidences
        logical(1)     :: column_is_predicted
        
        inf_predicted_state(:,:) = .FALSE.
        cell_confidence(:,:)     = 0
        col_confidence(:)        = 0
        
        n_predicted_columns = 0
        do column_i = 1, N_COLUMNS
            column_is_predicted = .FALSE.
            
            do segment_i = 1, columns(column_i) % dd_segment_count
                cell = columns(column_i) % dd_segment_destination(segment_i)
                if (cell /= -1) then
                    call get_segment_activity_level(                    &   
                        column                  = column_i,             &
                        segment                 = segment_i,            &
                        active_cells_all        = inf_active_state,     &
                        connected_synapses_only = .FALSE.,              &
                        ! out
                        n_active_synapses       = n_active_synapses)

                    if (n_active_synapses >= DD_ACTIVATION_THRESHOLD) then
                        dc = 1 ! TODO calc the duty cycles
                        cell_confidence(cell, column_i) = cell_confidence(cell, column_i) + dc
                        col_confidence(column_i) = col_confidence(column_i) + dc
                    
                        call get_segment_activity_level(                    &
                            column                  = column_i,             &
                            segment                 = segment_i,            &
                            active_cells_all        = inf_active_state,     &
                            connected_synapses_only = .TRUE.,               &
                            ! out
                            n_active_synapses       = n_active_synapses)

                        if (n_active_synapses >= DD_ACTIVATION_THRESHOLD) then
                            inf_predicted_state(cell, column_i) = .TRUE.
                            column_is_predicted = .TRUE.
                        end if
                    end if
                end if
            end do
            if (column_is_predicted) call inc(n_predicted_columns)
        end do

        ! normalize column and cell confidences
        sum_confidences = SUM(col_confidence(:))
        if (sum_confidences > 0.0) then
            col_confidence  = col_confidence  / sum_confidences
            cell_confidence = cell_confidence / sum_confidences
        end if
        
        ! are we predicting the required minimum number of columns?
        in_sequence = (n_predicted_columns >= (0.5 * avg_input_density))
        
    end subroutine infer_phase_2
!========================================================================
    subroutine infer_backtrack(                                         &
            active_columns,                                             &
            ! inout
            inf_predicted_state_prev,                                   &
            inf_active_state,                                           &
            ! out
            inf_predicted_state,                                        &
            col_confidence,                                             &
            cell_confidence)
        
        logical(1),   intent(in)  :: active_columns               (N_COLUMNS)
        logical(1),   intent(inout) :: inf_predicted_state_prev     (N_CELLS_PC, N_COLUMNS)
        logical(1),   intent(inout) :: inf_active_state             (N_CELLS_PC, N_COLUMNS)
        
        logical(1),   intent(out)    :: inf_predicted_state          (N_CELLS_PC, N_COLUMNS)
        real(REAL32), intent(out)   :: col_confidence               (N_COLUMNS)
        real(REAL32), intent(out)   :: cell_confidence              (N_CELLS_PC, N_COLUMNS)

    !This "backtracks" our inference state, trying to see if we can lock onto
    !the current set of inputs by assuming the sequence started up to N steps
    !ago on start cells.

    !@param activeColumns The list of active column indices

    !This will adjust @ref infActiveState['t'] if it does manage to lock on to a
    !sequence that started earlier. It will also compute infPredictedState['t']
    !based on the possibly updated @ref infActiveState['t'], so there is no need to
    !call inferPhase2() after calling inferBacktrack().
    !
    !This looks at:
    !    - @ref infActiveState['t']
    !
    !This updates/modifies:
    !    - @ref infActiveState['t']
    !    - @ref infPredictedState['t']
    !    - @ref colConfidence['t']
    !    - @ref cellConfidence['t']
    !
    !How it works:
    !-------------------------------------------------------------------
    !This method gets called from updateInferenceState when we detect either of
    !the following two conditions:
    !
    !-# The current bottom-up input had too many un-expected columns
    !-# We fail to generate a sufficient number of predicted columns for the
    !   next time step.
    !
    !Either of these two conditions indicate that we have fallen out of a
    !learned sequence.
    !
    !Rather than simply "giving up" and bursting on the unexpected input
    !columns, a better approach is to see if perhaps we are in a sequence that
    !started a few steps ago. The real world analogy is that you are driving
    !along and suddenly hit a dead-end, you will typically go back a few turns
    !ago and pick up again from a familiar intersection.
    !
    !This back-tracking goes hand in hand with our learning methodology, which
    !always tries to learn again from start cells after it loses context. This
    !results in a network that has learned multiple, overlapping paths through
    !the input data, each starting at different points. The lower the global
    !decay and the more repeatability in the data, the longer each of these
    !paths will end up being.
    !
    !The goal of this function is to find out which starting point in the past
    !leads to the current input with the most context as possible. This gives us
    !the best chance of predicting accurately going forward. Consider the
    !following example, where you have learned the following sub-sequences which
    !have the given frequencies:
    !
    !                  ? - Q - C - D - E      10X      seq 0
    !                  ? - B - C - D - F      1X       seq 1
    !                  ? - B - C - H - I      2X       seq 2
    !                  ? - B - C - D - F      3X       seq 3
    !          ? - Z - A - B - C - D - J      2X       seq 4
    !          ? - Z - A - B - C - H - I      1X       seq 5
    !          ? - Y - A - B - C - D - F      3X       seq 6
    !
    !        ----------------------------------------
    !      W - X - Z - A - B - C - D          <= input history
    !!                              ^
    !                              current time step
    !
    !Suppose, in the current time step, the input pattern is D and you have not
    !predicted D, so you need to backtrack. Suppose we can backtrack up to 6
    !steps in the past, which path should we choose? From the table above, we can
    !see that the correct answer is to assume we are in seq 4. How do we
    !implement the backtrack to give us this right answer? The current
    !implementation takes the following approach:
    !
    !-# Start from the farthest point in the past.
    !-# For each starting point S, calculate the confidence of the current
    !   input, conf(startingPoint=S), assuming we followed that sequence.
    !   Note that we must have learned at least one sequence that starts at
    !   point S.
    !-# If conf(startingPoint=S) is significantly different from
    !   conf(startingPoint=S-1), then choose S-1 as the starting point.
    !
    !The assumption here is that starting point S-1 is the starting point of
    !a learned sub-sequence that includes the current input in it's path and
    !that started the longest ago. It thus has the most context and will be
    !the best predictor going forward.
    !
    !From the statistics in the above table, we can compute what the confidences
    !will be for each possible starting point:
    !
    !    startingPoint           confidence of D
    !    -----------------------------------------
    !    B (t-2)               4/6  = 0.667   (seq 1,3)/(seq 1,2,3)
    !    Z (t-4)               2/3  = 0.667   (seq 4)/(seq 4,5)
    !
    !First of all, we do not compute any confidences at starting points t-1, t-3,
    !t-5, t-6 because there are no learned sequences that start at those points.
    !
    !Notice here that Z is the starting point of the longest sub-sequence leading
    !up to the current input. Event though starting at t-2 and starting at t-4
    !give the same confidence value, we choose the sequence starting at t-4
    !because it gives the most context, and it mirrors the way that learning
    !extends sequences.

        ! local variables
        logical(1)     :: inf_active_state_backup           (N_CELLS_PC, N_COLUMNS)
        logical(1)     :: inf_predicted_state_backup        (N_CELLS_PC, N_COLUMNS)
        integer(INT32) :: start_offset, current_time_step_offset, n_prev_patterns
        integer(INT32) :: offset
        real(REAL32)   :: total_confidence, cand_confidence
        integer(INT32) :: cand_start_offset
        
        logical(1)     :: inf_active_state_candidate        (N_CELLS_PC, N_COLUMNS)
        logical(1)     :: inf_predicted_state_candidate     (N_CELLS_PC, N_COLUMNS)
        real(REAL32)   :: cell_confidence_candidate         (N_CELLS_PC, N_COLUMNS)
        real(REAL32)   :: col_confidence_candidate          (N_COLUMNS)

        logical(1)     :: historic_active_columns           (N_COLUMNS)
        logical(1)     :: in_sequence
        logical(1)     :: use_start_cells
        logical(1)     :: bad_patterns                      (0:HISTORY_MAX_INF-1)

        
        n_prev_patterns = prev_inf_patterns_get_size()
        if (n_prev_patterns <= 0) then
            return !multiple procedure exits: ugly
        end if
        
        current_time_step_offset   = n_prev_patterns - 1
        
#       if _DEBUG
            if (.TRUE.) write(*, '("INFO:TP:infer_backtrack:n_prev_patterns"(1X,I0)"; current_time_step_offset"(1X,I0))'), n_prev_patterns, current_time_step_offset
            !read *
#       endif
        
        inf_active_state_backup    = inf_active_state
        inf_predicted_state_backup = inf_predicted_state_prev
        
    ! Let's go back in time and replay the recent inputs from start cells and
    !  see if we can lock onto this current set of inputs that way.
    !
    ! Start the farthest back and work our way forward. For each starting point,
    !  See if firing on start cells at that point would predict the current
    !  input as well as generate sufficient predictions for the next time step.
    !
    ! We want to pick the point closest to the current time step that gives us
    ! the relevant confidence. Think of this example, where we are at D and need
    ! to
    !   A - B - C - D
    ! decide if we should backtrack to C, B, or A. Suppose B-C-D is a high order
    ! sequence and A is unrelated to it. If we backtrack to B would we get a
    ! certain confidence of D, but if went went farther back, to A, the
    ! confidence wouldn't change, since A has no impact on the B-C-D series.
    !
    ! So, our strategy will be to pick the "B" point, since choosing the A point
    !  does not impact our confidences going forward at all.
        
   !    in_sequence = .FALSE.
        cand_confidence   = -1
        cand_start_offset = -1
            
        loop_a: do start_offset = 0, n_prev_patterns
            if ((start_offset == current_time_step_offset).AND.(cand_confidence /= -1)) then
                exit loop_a
            end if
            in_sequence = .FALSE.

            loop_b: do offset = start_offset, n_prev_patterns
                ! If we are about to set the active columns for the current time step
                ! based on what we predicted, capture and save the total confidence of
                ! predicting the current input
                if (offset == current_time_step_offset) then
                   total_confidence = SUM(col_confidence, MASK=active_columns)
                end if
                
                ! compute active_state given bottom up and predicted_state_prev
                inf_predicted_state_prev = inf_predicted_state

                historic_active_columns = prev_inf_patterns_get_data(offset)
                use_start_cells = (offset == start_offset)
                
                call infer_phase_1(                                             &
                    active_columns              = historic_active_columns,      &
                    inf_predicted_state_prev    = inf_predicted_state_prev,     &
                    use_start_cells             = use_start_cells,              &
                    ! out
                    inf_active_state            = inf_active_state,             &
                    in_sequence                 = in_sequence)

                if (NOT(in_sequence)) then
                    exit loop_b
                end if

#               if _DEBUG
                    if (.TRUE.) write(*, '("INFO:TP:infer_backtrack: before infer phase 2!")')
                    !read *
#               endif
                ! compute predicted_state given active_state
                call infer_phase_2(                                             &
                    inf_active_state        = inf_active_state,                 &
                    avg_input_density       = avg_input_density,                &
                    ! out
                    inf_predicted_state     = inf_predicted_state,              &
                    col_confidence          = col_confidence,                   &
                    cell_confidence         = cell_confidence,                  &
                    in_sequence             = in_sequence)

                if (NOT(in_sequence)) then
                    exit loop_b
                end if
            end do loop_b
            
            ! if starting from start_offset got lost along the way, mark it as an
            ! invalid start point
            if (NOT(in_sequence)) then
                bad_patterns(start_offset) = .TRUE.
            else
#               if _DEBUG
                    if (.TRUE.) write(*, '("INFO:TP:infer_backtrack:in sequence!")')
                    !read *
#               endif
                
                ! if we got to here, start_offset is a condidate starting point.
                ! save this state as a condidate state. It will become the chosen state if
                ! we detect a change in confidences starting at a later start_offset

                cand_confidence   = total_confidence
                cand_start_offset = start_offset

                if (cand_start_offset == current_time_step_offset) then ! no more to try
                    exit loop_a
                end if

                inf_active_state_candidate    = inf_active_state
                inf_predicted_state_candidate = inf_predicted_state
                cell_confidence_candidate     = cell_confidence
                col_confidence_candidate      = col_confidence
                exit loop_a
            end if
        end do loop_a

        ! if we failed to lock on at any starting point, fall back to the original
        ! active state that we had on entry
        if (cand_start_offset == -1) then
            inf_active_state = inf_active_state_backup
            call infer_phase_2(                                         &
                inf_active_state        = inf_active_state,             &
                avg_input_density       = avg_input_density,            &
                ! out
                inf_predicted_state     = inf_predicted_state,          &
                col_confidence          = col_confidence,               &
                cell_confidence         = cell_confidence,              &
                in_sequence             = in_sequence)
            
        else
            ! install the candidate state, it it was not the last we evaluated
            if (cand_start_offset /= current_time_step_offset) then
                inf_active_state    = inf_active_state_candidate
                inf_predicted_state = inf_predicted_state_candidate
                cell_confidence     = cell_confidence
                col_confidence      = col_confidence
            end if
        end if

        ! remove any useless patterns at the head of the previous input pattern
        ! queue.
        do offset = 1, prev_inf_patterns_get_size()
            if (bad_patterns(offset)) then
                !TODO
                !call prev_inf_patterns_remove_last()
            end if
        end do
        
        ! restore the original predicted state
        inf_predicted_state_prev = inf_predicted_state_backup
    end subroutine infer_backtrack
!========================================================================
!========================================================================
    subroutine update_learning_state(                                   &
            active_columns,                                             &
            lrn_predicted_state_prev,                                   &
            lrn_active_state_prev,                                      &
            ! out
            lrn_predicted_state,                                        &
            lrn_active_state)

        logical(1), intent(in)  :: active_columns               (N_COLUMNS)
        logical(1), intent(inout)  :: lrn_predicted_state_prev     (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(inout)  :: lrn_active_state_prev        (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: lrn_predicted_state          (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: lrn_active_state             (N_CELLS_PC, N_COLUMNS)
        ! local variables
        logical(1)     :: in_sequence
        integer(INT32) :: column_i, back_steps
        

        ! update our learning input history
        call prev_lrn_patterns_add(active_columns)
        
        ! process queued up segment updates, now that we have bottom-up, we
        ! can update the permanences on the cells that we predicted to to turn on
        ! and did receive bottom-up
        do column_i = 1, N_COLUMNS
            if (active_columns(column_i)) then
                call process_segment_updates(                       &
                    column              = column_i,                 &
                    column_active       = .TRUE.,                   &
                    lrn_iteration_idx   = lrn_iteration_idx)
            end if
        end do

        ! Decrement the PAM counter if it is running 
        if (pam_counter > 0) pam_counter = pam_counter - 1
        ! increment our learned sequence length
        call inc(learned_seq_length)

        ! phase 1 - turn on the predicted cell in the column that received
        ! bottom-up. If there was no predicted cell, pick one to learn on
        if (NOT(reset_called)) then
            call learn_phase_1(                                             &
                active_columns              = active_columns,               &
                lrn_predicted_state_prev    = lrn_predicted_state_prev,     &
                lrn_active_state_prev       = lrn_active_state_prev,        &
                read_only                   = .FALSE.,                      &
                ! out
                lrn_active_state            = lrn_active_state,             &
                in_sequence                 = in_sequence)

            ! reset our PAM counter if we are in sequence
            if (in_sequence) pam_counter = pam_length
        end if

        ! start over on start cells if any of the following occurs:
        ! 1] a reset was just called
        ! 2] we have been too long out of sequence (the pam_counter has expired)
        ! 3] we have reached maximum allowed sequence length
        !
        ! Note that, unless we are following a reset, we als just learned or
        ! re-enforced connections to the current set of active columns because
        ! this input is still a valid prediction to learn
        !
        ! It is especially helpful to learn the connections to this input when
        ! you have a max_seq_length constraint in place. Otherwise, you will have
        ! no continuity at all between sub-sequences of length max_seq_length.
        
        if ((reset_called).OR.                                      &
            (pam_counter == 0).OR.(                                 &
                (MAX_SEQUENCE_LENGTH /= 0).AND.                     &
                (learned_seq_length >= MAX_SEQUENCE_LENGTH))) then
            
            ! update average learned sequence length - this is a diagnostic statistic
            if (pam_counter == 0) then
                call update_avg_learned_seq_length(learned_seq_length - pam_length)
            else
                call update_avg_learned_seq_length(learned_seq_length)
            end if
            
            ! backtrack to an earlier starting point, if we find one
            back_steps = 0
            if (NOT(reset_called)) then
                call learn_backtrack(                                       &
                    back_steps                  = back_steps,               &
                    lrn_predicted_state         = lrn_predicted_state,      &
                    lrn_predicted_state_prev    = lrn_predicted_state_prev, &
                    lrn_active_state            = lrn_active_state,         &
                    lrn_active_state_prev       = lrn_active_state_prev)
            end if

            ! start over in the current time step if reset was called, we we could not
            ! backtrack
            if ((reset_called).OR.(back_steps == 0)) then
                lrn_active_state(1, :) = active_columns(:)
                lrn_predicted_state = .FALSE.
                ! remove any old input history patterns
                call prev_lrn_patterns_clear()
            end if

            ! reset PAM counter
            pam_counter        = pam_length
            learned_seq_length = back_steps

            ! clear out any old segment updates from prior sequences
            do column_i = 1, N_COLUMNS
                call clear_segment_updates(column_i)
            end do
        end if
        
        ! Phase2: compute new predicted state. When computing predictions for 
        ! phase 2, we predict at most one cell per columns (the one with the best
        ! matching segment
        call learn_phase_2(                                             &
            read_only               = .FALSE.,                          &
            lrn_active_state        = lrn_active_state,                 &
            lrn_active_state_prev   = lrn_active_state_prev,            &
            ! out
            lrn_predicted_state     = lrn_predicted_state)
        
    end subroutine update_learning_state
!========================================================================
    subroutine update_avg_learned_seq_length(prev_seq_length)
        integer(INT32), intent(in) :: prev_seq_length
        real(REAL32) :: alpha
        if (lrn_iteration_idx < 100) then
            alpha = 0.5
        else
            alpha = 0.1
        end if
        avg_learned_seq_length = (((1.0 - alpha) * avg_learned_seq_length) + (alpha * prev_seq_length))
        
    end subroutine update_avg_learned_seq_length
!========================================================================
    subroutine learn_phase_1(                                           &
            active_columns,                                             &
            lrn_predicted_state_prev,                                   &
            lrn_active_state_prev,                                      &
            read_only,                                                  &
            ! out
            lrn_active_state,                                           &
            in_sequence)

    ! read_only: true if being called from backtracking logic. THis tells us not
    ! to increment any segment duty cycles or queue up any updates
    ! in_sequence: is set to true if the current input was sufficiently predicted, 
    ! OR if we started over on start_cells. False indicates that the current 
    ! input was NOT predicted, well enough to consider it as "in_sequence" .
    
        logical(1), intent(in)  :: active_columns               (N_COLUMNS)
        logical(1), intent(in)  :: lrn_predicted_state_prev     (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(in)  :: lrn_active_state_prev        (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(in)  :: read_only
        logical(1), intent(out) :: lrn_active_state             (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: in_sequence
        ! local variables
        integer(INT32) :: column_i, num_unpredicted_columns, unused, best_segment
        integer(INT8)  :: cell_i, predicted_cell, best_cell
        type(seg_update_t) :: seg_update
        
        ! save previous active state and start out on a clean slate
        lrn_active_state = .FALSE.

        ! for each column, turn on the predicted cell. There will always be at most
        ! one predicted cell per column
        num_unpredicted_columns = 0
        column_loop: do column_i = 1, N_COLUMNS
            if (active_columns(column_i)) then
                
                predicted_cell = -1
                do cell_i = 1, N_CELLS_PC
                    if (lrn_predicted_state_prev(cell_i, column_i)) then
                        predicted_cell = cell_i
                    end if
                end do

                ! if we have a predicted cell, turn it on. The segment's pos_activation
                ! count will have already been incremented by process_segment_updates

                if (predicted_cell /= -1) then
                    lrn_active_state(predicted_cell, column_i) = .TRUE.
                    cycle column_loop
                end if
                
                call inc(num_unpredicted_columns)
                
                if (read_only) then
                    cycle column_loop
                end if

                ! if no predicted cell, pick the closest matching one to reinforce, or
                ! if none exists, create a new segment on a cell in that column
                call get_best_matching_cell(                            &
                    column              = column_i,                     &
                    active_cells_all    = lrn_active_state_prev,        &
                    threshold           = MIN_DD_ACTIVATION_THRESHOLD,  &
                    ! out
                    best_cell           = best_cell,                    &
                    best_segment        = best_segment,                 &
                    best_activity       = unused)
                
                if (best_segment /= -1) then
                    if (columns(column_i) % sequence_segments(best_segment)) then
                        lrn_active_state(best_cell, column_i) = .TRUE.

                        call get_segment_active_synapses(                   &
                            column              = column_i,                 &
                            cell                = best_cell,                &
                            segment             = best_segment,             &
                            active_cells_all    = lrn_active_state_prev,    &
                            new_synapses        = .TRUE.,                   &
                            ! out
                            seg_update          = seg_update)
                    
                        call inc(columns(column_i) % dd_synapse_total_activation(best_segment))
                        call adapt_segment(                                 &
                            column              = column_i,                 &
                            seg_update          = seg_update,               &
                            lrn_iteration_idx   = lrn_iteration_idx,        &
                            permanence_inc      = DD_PERMANENCE_INC,        &
                            permanence_dec      = DD_PERMANENCE_DEC)
                    end if
                else 
                    ! if no close matching segment exists, create a new segment
                    ! choose a cell in this column to add a new segment to
                    call get_cell_for_new_segment(cell=best_cell, column=column_i)
                    lrn_active_state(best_cell, column_i) = .TRUE.

                    call get_segment_active_synapses(                   &
                        column              = column_i,                 &
                        cell                = best_cell,                &
                        segment             = -1,                       &
                        active_cells_all    = lrn_active_state_prev,    &
                        new_synapses        = .TRUE.,                   &
                        ! out
                        seg_update          = seg_update)

                    seg_update % sequence_segment = .TRUE.
                    call adapt_segment(                                 &
                        column              = column_i,                 &
                        seg_update          = seg_update,               &
                        lrn_iteration_idx   = lrn_iteration_idx,        &
                        permanence_inc      = DD_PERMANENCE_INC,        &
                        permanence_dec      = DD_PERMANENCE_DEC)
                end if
            end if
        end do column_loop

        ! determine if we are out of sequences or not
        in_sequence = (num_unpredicted_columns < (count_active(active_columns) / 2))
        
    end subroutine learn_phase_1
!========================================================================
    subroutine learn_phase_2(                                           &
            read_only,                                                  &
            lrn_active_state,                                           &
            lrn_active_state_prev,                                      &
            ! out
            lrn_predicted_state)

    ! this computes the lrn_predicted_state and queues up any segments that
    ! became active (and the list of active synapses for each segment) into
    ! the segment_updates queueu

        logical(1), intent(in)  :: read_only
        logical(1), intent(in)  :: lrn_active_state             (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(in)  :: lrn_active_state_prev        (N_CELLS_PC, N_COLUMNS)
        logical(1), intent(out) :: lrn_predicted_state          (N_CELLS_PC, N_COLUMNS)
        
        ! local variables
        integer(INT32) :: column_i, total_activation, best_activity, best_segment
        integer(INT8)  :: best_cell
        logical(1)     :: new_synapses
        type(seg_update_t) :: seg_update
        
        ! clear out predicted state to start with
        lrn_predicted_state = .FALSE.

        ! compute new predicte state. When computing predictions for
        ! phase 2, we predict at most one cell per column (the one with the best
        ! matching segment).
        column_loop: do column_i = 1, N_COLUMNS

            ! is there a cell predicted to turn on in this column?
            call get_best_matching_cell(                                &
                column              = column_i,                         &
                active_cells_all    = lrn_active_state,                 &
                threshold           = DD_ACTIVATION_THRESHOLD,          &
                ! out
                best_cell           = best_cell,                        &
                best_segment        = best_segment,                     &
                best_activity       = best_activity)

            if (best_cell == -1) then
                cycle column_loop
            end if

            ! turn on the predicted state for the best matching cell and queue
            ! the pertinent segment up for an update, which will get processed if
            ! the cell receives bottom up in the future.
            lrn_predicted_state(best_cell, column_i) = .TRUE.

            if (read_only) then
                cycle column_loop
            end if

            ! queue up this segment for updating
            new_synapses = (best_activity < NEW_NUMBER_SYNAPSES)
            call get_segment_active_synapses(                   &
                column              = column_i,                 &
                cell                = best_cell,                &
                segment             = best_segment,             &
                active_cells_all    = lrn_active_state,         &
                new_synapses        = new_synapses,             &
                ! out
                seg_update          = seg_update)

            call inc(total_activation)
            call add_to_segment_updates(column_i, seg_update)
            
            if (DO_POOLING) then
                ! creates a new pooling segment if no best matching found
                ! sum(all synapses) >= min_threshold, "weak" activations
                call get_best_matching_segment(                     &
                    column              = column_i,                 &
                    cell                = best_cell,                &
                    active_cells_all    = lrn_active_state_prev,    &
                    ! out
                    best_segment        = best_segment)
                
                call get_segment_active_synapses(                   &
                    column              = column_i,                 &
                    cell                = best_cell,                &
                    segment             = best_segment,             &
                    active_cells_all    = lrn_active_state_prev,    &
                    new_synapses        = .TRUE.,             &
                    ! out
                    seg_update          = seg_update)
                
                call add_to_segment_updates(column_i, seg_update)
            end if
        end do column_loop

    end subroutine learn_phase_2
!========================================================================
    subroutine learn_backtrack(                                         &
            ! out
            lrn_predicted_state,                                        &
            lrn_predicted_state_prev,                                   &
            lrn_active_state_prev,                                      &
            lrn_active_state,                                           &
            back_steps)
    
        logical(1),     intent(out) :: lrn_predicted_state_prev   (N_CELLS_PC, N_COLUMNS)
        logical(1),     intent(out) :: lrn_active_state_prev      (N_CELLS_PC, N_COLUMNS)
        logical(1),     intent(out)   :: lrn_active_state           (N_CELLS_PC, N_COLUMNS)
        logical(1),     intent(out)   :: lrn_predicted_state        (N_CELLS_PC, N_COLUMNS)
        integer(INT32), intent(out)   :: back_steps
        
        ! This "backtracks" our learning state, trying to see if we can lock onto
        ! the current set of inputs by assuming the sequence started up to N steps
        ! ago on start cells.

        ! This will adjust @ref lrnActiveState['t'] if it does manage to lock on to a
        ! sequence that started earlier.

        !@returns          >0 if we managed to lock on to a sequence that started
        !                 earlier. The value returned is how many steps in the
        !                 past we locked on.
        !                 If 0 is returned, the caller needs to change active
        !                 state to start on start cells.

        !How it works:
        !-------------------------------------------------------------------
        !This method gets called from updateLearningState when we detect either of
        !the following two conditions:
    
        !-# Our PAM counter (@ref pamCounter) expired
        !-# We reached the max allowed learned sequence length

        !Either of these two conditions indicate that we want to start over on start
        !cells.
        !
        !Rather than start over on start cells on the current input, we can
        !accelerate learning by backtracking a few steps ago and seeing if perhaps
        !a sequence we already at least partially know already started.

        !This updates/modifies:
        !    - @ref lrnActiveState['t']

        !This trashes:
        !    - @ref lrnActiveState['t-1']
        !    - @ref lrnPredictedState['t']
        !    - @ref lrnPredictedState['t-1']

    
        ! local variables
        logical(1)     :: in_sequence
        integer(INT32) :: num_prev_patterns, i, start_offset
        logical(1)     :: bad_patterns          (0:HISTORY_MAX_LRN-1)
    
        ! How much input history have we accumulated?
        ! The current input is always at the end of self._prevInfPatterns (at
        ! index -1), and is not a valid startingOffset to evaluate.
        num_prev_patterns = prev_lrn_patterns_get_size()
        if (num_prev_patterns <= 0) then
            back_steps = 0
            return
        end if

        ! We will record which previous input patterns did not generate predictions
        ! up to the current time step and remove all the ones at the head of the
        ! input history queue so that we don't waste time evaluating them again at
        ! a later time step.
        bad_patterns(:) = .FALSE.

        ! Let's go back in time and replay the recent inputs from start cells and
        ! see if we can lock onto this current set of inputs that way.
        !
        ! Start the farthest back and work our way forward. For each starting point,
        ! See if firing on start cells at that point would predict the current
        ! input.
        !
        ! We want to pick the point farthest in the past that has continuity
        ! up to the current time step
        start_offset = -1
        
        loop_a: do i = 0, num_prev_patterns
            ! Can we backtrack from startOffset?
            call learn_backtrack_from(                                  &
                start_offset                = i,                        &
                read_only                   = .TRUE.,                   &
                ! inout
                lrn_active_state            = lrn_active_state,         &
                lrn_predicted_state         = lrn_predicted_state,      &
                ! out
                lrn_active_state_prev       = lrn_active_state_prev,    &
                lrn_predicted_state_prev    = lrn_predicted_state_prev, &
                in_sequence                 = in_sequence)

            ! Done playing through the sequence from starting point startOffset
            ! Break out as soon as we find a good path
            if (in_sequence) then
                start_offset = i
                exit loop_a
            else
                ! Take this bad starting point out of our input history so we don't
                ! try it again later.
                bad_patterns(i) = .TRUE.
            end if
        end do loop_a
        
        ! If we failed to lock on at any starting point, return failure. The caller
        ! will start over again on start cells
        if (start_offset /= -1) then
            ! Nothing in our input history was a valid starting point, so get rid
            ! of it so we don't try any of them again at a later iteration
            call prev_lrn_patterns_clear()
            back_steps = 0
        else
            ! We did find a valid starting point in the past. Now, we need to
            ! re-enforce all segments that became active when following this path.
 
            call learn_backtrack_from(                                  &
                start_offset                = start_offset,             &
                read_only                   = .FALSE.,                  &
                ! inout
                lrn_active_state            = lrn_active_state,         &
                lrn_predicted_state         = lrn_predicted_state,      &
                ! out
                lrn_active_state_prev       = lrn_active_state_prev,    &
                lrn_predicted_state_prev    = lrn_predicted_state_prev, &
                in_sequence                 = in_sequence)

            ! Remove any useless patterns at the head of the input pattern history queue.
            loop_b: do i = 0, num_prev_patterns
                if (bad_patterns(i).OR.(i <= start_offset)) then
                    !TODO
                    !call prev_lrn_patterns_remove_last()
                else
                    exit loop_b
                end if
            end do loop_b
            back_steps = num_prev_patterns - start_offset
        end if
    end subroutine learn_backtrack
!========================================================================
    subroutine learn_backtrack_from(                                    &
            start_offset,                                               &
            read_only,                                                  &
            ! inout
            lrn_active_state,                                           &
            lrn_predicted_state,                                        &
            ! out
            lrn_predicted_state_prev,                                   &
            lrn_active_state_prev,                                      &
            in_sequence)
    
    ! A utility method called from learnBacktrack. This will backtrack
    ! starting from the given startOffset in our prevLrnPatterns queue.

    ! It returns True if the backtrack was successful and we managed to get
    ! predictions all the way up to the current time step.

    ! If readOnly, then no segments are updated or modified, otherwise, all
    ! segment updates that belong to the given path are applied.

    !@param startOffset Start offset within the prevLrnPatterns input history
    !@returns           True if we managed to lock on to a sequence that started
    !                   earlier.
    !                  If False, we lost predictions somewhere along the way
    !                   leading up to the current time.
 
        integer(INT32), intent(in)    :: start_offset
        logical(1),     intent(in)    :: read_only

        logical(1),     intent(inout)   :: lrn_active_state           (N_CELLS_PC, N_COLUMNS)
        logical(1),     intent(inout)   :: lrn_predicted_state        (N_CELLS_PC, N_COLUMNS)

        logical(1),     intent(out)   :: lrn_predicted_state_prev   (N_CELLS_PC, N_COLUMNS)
        logical(1),     intent(out)   :: lrn_active_state_prev      (N_CELLS_PC, N_COLUMNS)
        logical(1),     intent(out)   :: in_sequence
        ! local variables
        integer(INT32) :: num_prev_patterns, current_time_steps_offset, offset, column_i
        logical(1)     :: input_columns     (N_COLUMNS)
        
        in_sequence = .FALSE.

        ! How much input history have we accumulated?
        ! The current input is always at the end of self._prevInfPatterns (at
        ! index -1), but it is also evaluated as a potential starting point by
        ! turning on it's start cells and seeing if it generates sufficient
        ! predictions going forward.
        num_prev_patterns = prev_lrn_patterns_get_size()

        ! This is an easy to use label for the current time step
        current_time_steps_offset = num_prev_patterns - 1

        ! Clear out any old segment updates. learnPhase2() adds to the segment
        ! updates if we're not readOnly
        if (NOT(read_only)) then
            do column_i = 1, N_COLUMNS
                call clear_segment_updates(column_i)
            end do
        end if
    
        ! Play through up to the current time step
        in_sequence = .TRUE.
        do offset = start_offset, num_prev_patterns

            ! Copy predicted and active states into t-1
            lrn_predicted_state_prev = lrn_predicted_state
            lrn_active_state_prev    = lrn_active_state

            ! Get the input pattern
            input_columns = prev_lrn_patterns_get_data(offset)

            ! Apply segment updates from the last set of predictions
            if (NOT(read_only)) then
                do column_i = 1, N_COLUMNS
                    if (input_columns(column_i)) then
                        call process_segment_updates(               &
                            column              = column_i,         &
                            column_active       = .TRUE.,           &
                            lrn_iteration_idx   = lrn_iteration_idx)
                    end if
                end do
            end if

            ! Phase 1:
            ! Compute activeState[t] given bottom-up and predictedState[t-1]
            if (offset == start_offset) then
                lrn_active_state = .FALSE.
                lrn_active_state(1, :) = input_columns(:)
                in_sequence = .TRUE.
            else
                ! Uses lrnActiveState['t-1'] and lrnPredictedState['t-1']
                ! computes lrnActiveState['t']
                call learn_phase_1(                                             &
                    active_columns              = input_columns,                &
                    lrn_predicted_state_prev    = lrn_predicted_state_prev,     &
                    lrn_active_state_prev       = lrn_active_state_prev,        &
                    read_only                   = read_only,                    &
                    ! out
                    lrn_active_state            = lrn_active_state,             &
                    in_sequence                 = in_sequence)
            end if
      
            ! Break out immediately if we fell out of sequence or reached the current
            ! time step
            if (NOT(in_sequence).OR.(offset == current_time_steps_offset)) then
                exit
            end if
            
            ! Phase 2:
            ! Computes predictedState['t'] given activeState['t'] and also queues
            ! up active segments into self.segmentUpdates, unless this is readOnly
            call learn_phase_2(                                             &
                read_only               = read_only,                        &
                lrn_active_state        = lrn_active_state,                 &
                lrn_active_state_prev   = lrn_active_state_prev,            &
                ! out
                lrn_predicted_state     = lrn_predicted_state)
        end do
        ! inSequence is whether or not this starting point was valid
    end subroutine
!========================================================================
!========================================================================
    subroutine apply_decay()
        !TODO
    end subroutine apply_decay
!========================================================================
    subroutine get_segment_activity_level(                              &
            column,                                                     &
            segment,                                                    &
            active_cells_all,                                           &
            connected_synapses_only,                                    &
            !out
            n_active_synapses)
            
        integer(INT32), intent(in)  :: column
        integer(INT32), intent(in)  :: segment
        logical(1),     intent(in)  :: active_cells_all     (N_CELLS)
        logical(1),     intent(in)  :: connected_synapses_only
        integer(INT32), intent(out) :: n_active_synapses
        ! local variables
        integer(INT32) :: synapse_i, global_cell_id

        n_active_synapses = 0
        
        if (connected_synapses_only) then
            do synapse_i = 1, N_DD_SYNAPSES
                global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment)
                if (global_cell_id == -1) exit
                if (is_DD_connected(column, segment, synapse_i, DD_PERMANENCE_THRESHOLD)) then
                    if (active_cells_all(global_cell_id)) then
                       call inc(n_active_synapses)
                    end if
                end if
            end do
        else
            do synapse_i = 1, N_DD_SYNAPSES
                global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment)
                if (global_cell_id == -1) exit
                if (active_cells_all(global_cell_id)) then
                    call inc(n_active_synapses)
                end if
            end do
        end if
        
#       if _DEBUG
            if (.FALSE.) print '("INFO:TP:get_segment_activity_level: column"(1X,I0)"; segment"(1X,I0)"; n_active_synapses"(1X,I0),"; N_DD_SYNAPSES",(1X,I0))', column, segment, n_active_synapses, N_DD_SYNAPSES
#       endif
    end subroutine get_segment_activity_level
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
    subroutine get_best_matching_cell(                                  &
            column,                                                     &
            active_cells_all,                                           &
            threshold,                                                  &
            ! out
            best_cell,                                                  &
            best_segment,                                               &
            best_activity)
    
    ! Find weakly activated cell in column with at least minThreshold active synapses.
    ! threshold minimum number of synapses required

        integer(INT32), intent(in)  :: column
        logical(1),     intent(in)  :: active_cells_all         (N_CELLS)
        integer(INT32), intent(in)  :: threshold
        integer(INT8),  intent(out) :: best_cell
        integer(INT32), intent(out) :: best_segment
        integer(INT32), intent(out) :: best_activity
        ! local variables
        integer(INT32) :: segment_i, activity

        ! Collect all cells in column c that have at least minThreshold in the most
        ! activated segment
        best_activity = threshold
        best_segment  = -1

        do segment_i = 1, columns(column) % dd_segment_count
            call get_segment_activity_level(                            &
                column                      = column,                   &
                segment                     = segment_i,                &
                active_cells_all            = active_cells_all,         &
                connected_synapses_only     = .TRUE.,                   &
                !out
                n_active_synapses           = activity)

            if (activity > best_activity) then
                best_activity = activity
                best_segment  = segment_i
            end if
        end do

        if (best_segment == -1) then
            best_cell           = -1
            best_segment        = -1
            best_activity       = -1
        else
            best_cell           = columns(column) % dd_segment_destination(best_segment)
            !best_segment is already set
            !best_activity is aready set
        end if
    end subroutine get_best_matching_cell
!========================================================================
    subroutine get_best_matching_segment(                               &
        column,                                                         &
        cell,                                                           &
        active_cells_all,                                               &
        ! out
        best_segment)
    
    ! For the given cell, find the segment with the largest number of active
    ! synapses. This routine is aggressive in finding the best match. The
    ! permanence value of synapses is allowed to be below connectedPerm. The number
    ! of active synapses is allowed to be below activationThreshold, but must be
    ! above minThreshold. The routine returns the segment index. If no segments are
    ! found, then an index of -1 is returned.

        integer(INT32), intent(in)  :: column
        integer(INT8),  intent(in)  :: cell
        logical(1),     intent(in)  :: active_cells_all         (N_CELLS)
        integer(INT32), intent(out) :: best_segment
        ! local variables
        integer(INT32) :: activity, best_activity, segment_i
        
        best_activity = MIN_DD_ACTIVATION_THRESHOLD
        best_segment  = -1
        
        do segment_i = 1, columns(column) % dd_segment_count
            if (columns(column) % dd_segment_destination(segment_i) == cell) then
                call get_segment_activity_level(                        &
                    column                  = column,                   &
                    segment                 = segment_i,                &
                    active_cells_all        = active_cells_all,         &
                    connected_synapses_only = .FALSE.,                  &
                    n_active_synapses       = activity)
                
                if (activity > best_activity) then
                    best_activity = activity
                    best_segment  = segment_i
                end if
            end if
        end do
    end subroutine get_best_matching_segment
!========================================================================
    subroutine trim_segments(column)
        integer(INT32), intent(in) :: column
        !NOT USED
    end subroutine
!========================================================================
    subroutine get_segment_active_synapses(                             &
        column,                                                         &
        cell,                                                           &
        segment,                                                        &
        active_cells_all,                                               &
        new_synapses,                                                   &
        !out
        seg_update) 
    
    !Return a segmentUpdate data structure containing a list of proposed
    !changes to segment s. Let activeSynapses be the list of active synapses
    !where the originating cells have their activeState output = 1 at time step
    !t. (This list is empty if s is None since the segment doesn't exist.)
    !newSynapses is an optional argument that defaults to false. If newSynapses
    !is true, then newSynapseCount - len(activeSynapses) synapses are added to
    !activeSynapses. These synapses are randomly chosen from the set of cells
    !that have learnState = 1 at timeStep.

        integer(INT32), intent(in) :: column
        integer(INT8),  intent(in) :: cell
        integer(INT32), intent(in) :: segment
        logical(1),     intent(in) :: active_cells_all      (N_CELLS)
        logical(1),     intent(in) :: new_synapses
        type(seg_update_t), intent(out) :: seg_update
        ! local variables
        integer(INT32) :: synapse_i, global_cell_id
        
        seg_update % dd_synapse_active = .FALSE.
        seg_update % cell = cell

        if (segment /= -1) then
            do synapse_i = 1, N_DD_SYNAPSES
                global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment)
                if (global_cell_id == -1) exit
                if (active_cells_all(global_cell_id)) then
                    seg_update % dd_synapse_active(synapse_i) = .TRUE.
                end if
            end do
        end if
        
        if (new_synapses) then ! add a few more synapses
            !TODO
            !n_synapses_to_add = self.newSynapseCount - len(activeSynapses)

            ! Here we add *pairs* (colIdx, cellIdx) to activeSynapses
            !activeSynapses += self.chooseCellsToLearnFrom(c, i, s, nSynapsesToAdd, activeState)

            ! It's still possible that activeSynapses is empty, and this will
            ! be handled in addToSegmentUpdates

            ! NOTE: activeSynapses contains a mixture of integers and pairs of integers
            ! - integers are indices of synapses already existing on the segment,
            !   that we will need to update.
            ! - pairs represent source (colIdx, cellIdx) of new synapses to create on
            !   the segment
        end if
        end subroutine
!========================================================================
    subroutine get_cell_for_new_segment(                                &
            column,                                                     &
            ! out
            cell)
    
        !Return the index of a cell in this column which is a good candidate
        !for adding a new segment.
    
        integer(INT32), intent(in)  :: column
        integer(INT8),  intent(out) :: cell
        
        !if we have fixed size resources in effect, we insure that we pick a
        !cell which does not already have the max number of allowed segments. If
        !none exists, we choose the least used segment in the column to re-allocate.

        ! NOTE: It is important NOT to always pick the cell with the fewest number
        ! of segments. The reason is that if we always do that, we are more likely
        ! to run into situations where we choose the same set of cell indices to
        ! represent an 'A' in both context 1 and context 2. This is because the
        ! cell indices we choose in each column of a pattern will advance in
        ! lockstep (i.e. we pick cell indices of 1, then cell indices of 2, etc.).

        if (N_CELLS_PC > 1) then
            ! Don't ever choose the start cell (cell # 1) in each column
            cell = rand_int32(2, N_CELLS_PC)
            
#           if _DEBUG
                if (.FALSE.) print '("INFO:TP:get_cell_for_new_segment: random cell =",(1X,I0))', cell
#           endif
        else
            cell = 0
        end if
    end subroutine
!========================================================================
    subroutine duty_cycles(                                             &
            column,                                                     &
            segment,                                                    &
            active,                                                     &
            read_only)
    
        integer(INT32), intent(in) :: column
        integer(INT32), intent(in) :: segment
        logical(1),     intent(in) :: active
        logical(1),     intent(in) :: read_only
        
        if (NOT(read_only)) then
            
        end if
    end subroutine
!========================================================================
    subroutine update_stats_infer_end(                                  &
!            internal_stats,                                             &
            active_columns,                                             &
!            predicted_stats,                                            &
            col_confidence)
        
!        internal_stats
        logical(1), intent(in)  :: active_columns       (N_COLUMNS)
        real(REAL32), intent(out) :: col_confidence     (N_COLUMNS)
        col_confidence = 0
        !TODO
    end subroutine update_stats_infer_end
!========================================================================
    subroutine inc(i)
        integer(INT32), intent(inout) :: i
        i = i + 1
    end subroutine
!========================================================================
!========================================================================
    subroutine clear_segment_updates(                                   &
            column)
        integer(INT32), intent(in) :: column
        columns(column) % seg_updates(:) % cell = -1
    end subroutine
!========================================================================
    subroutine add_to_segment_updates(                                  &
            column,                                                     &
            segment_update)
    
        integer(INT32),     intent(in) :: column
        type(seg_update_t), intent(in) :: segment_update
        ! local variables
        integer(INT32) :: i
    
#       if _DEBUG
            if (segment_update % cell == -1) print '("ERROR:TP:add_to_segment_updates: column"(1X,I0)": incorrect cell")', column
#       endif
        do i = 1, N_SEGMENT_UPDATES_MAX
            if (columns(column) % seg_updates(i) % cell /= -1) then
                columns(column) % seg_updates(i) = segment_update
            end if
            return
        end do
        print '("WARNING:TP:add_to_segment_updates: column"(1X,I0)": no empty segment slot present, ignoring the update")', column
    end subroutine
!========================================================================
    subroutine process_segment_updates(                                 &
            column,                                                     &
            column_active,                                              &
            lrn_iteration_idx)

!    Go through the list of accumulated segment updates and process them
!    as follows:

!    if the segment update is too old, remove the update,
!    else if the cell received bottom-up, update its permanences
!    else if it's still being predicted, leave it in the queue
!    else remove it.

        integer(INT32), intent(in) :: column
        logical(1),     intent(in) :: column_active
        integer(INT32), intent(in) :: lrn_iteration_idx
        ! local variables
        integer(INT32) :: i
        
        if (column_active) then
            do i = 1, N_SEGMENT_UPDATES_MAX
                if (columns(column) % seg_updates(i) % cell /= -1) then
                    ! If this segment has expired. Ignore this update (and hence remove it from list)
                    if ((lrn_iteration_idx - columns(column) % seg_updates(i) % create_date) > SEG_UPDATE_VALID_DURATION) then
                        ! the update has expired, remove it
                         columns(column) % seg_updates(i) % cell = -1
                    else
                        call adapt_segment(                                 &
                            column              = column,                   &
                            seg_update          = columns(column) % seg_updates(i), &
                            lrn_iteration_idx   = lrn_iteration_idx,        &
                            permanence_inc      = DD_PERMANENCE_INC,        &
                            permanence_dec      = DD_PERMANENCE_DEC)
                    end if
                end if
            end do
        end if
    end subroutine process_segment_updates
!========================================================================
    subroutine adapt_segment(                                           &
            column,                                                     &
            seg_update,                                                 &
            lrn_iteration_idx,                                          &
            permanence_inc,                                             &
            permanence_dec)
        
        !This function applies segment update information to a segment in a cell.
        !Synapses on the active list get their permanence counts incremented by
        !permanenceInc. All other synapses get their permanence counts decremented
        !by permanenceDec.

        !We also increment the positiveActivations count of the segment.

        integer(INT32),     intent(in) :: column
        type(seg_update_t), intent(in) :: seg_update
        integer(INT32),     intent(in) :: lrn_iteration_idx
        real(REAL32),       intent(in) :: permanence_inc
        real(REAL32),       intent(in) :: permanence_dec
        ! local variables
        real(REAL32)   :: old_permanence, new_permanence
        integer(INT32) :: synapse_i, global_cell_id, segment
        
        segment = seg_update % segment
        
        if (segment /= -1) then ! Modify an existing segment

            ! Mark it as recently useful
            columns(column) % dd_synapse_last_active_iteration(segment) = lrn_iteration_idx

            ! Update frequency and positiveActivations
            call inc(columns(column) % dd_synapse_positive_activations(segment))
            call duty_cycles(               &
                column      = column,       &
                segment     = segment,      &
                active      = .TRUE.,       &
                read_only   = .FALSE.)

            do synapse_i = 1, N_DD_SYNAPSES
                global_cell_id = columns(column) % dd_synapse_origin(synapse_i, segment)
                if (global_cell_id /= -1) then
                    old_permanence  = columns(column) % dd_synapse_permanence(synapse_i, segment)
                    if (seg_update % dd_synapse_active(synapse_i)) then
                        new_permanence = MIN(1.0, old_permanence + permanence_inc)
                    else
                        new_permanence = MAX(0.0, old_permanence - permanence_dec)
                    end if
                    columns(column) % dd_synapse_permanence(synapse_i, segment) = new_permanence
                else
                    exit
                end if
            end do
            
            ! Finally, create new synapses if needed
            if (seg_update % n_new_synapes > 0) then
                ! If we have fixed resources, get rid of some old syns if necessary
                ! if ((self.maxSynapsesPerSegment > 0).AND.(len(synsToAdd) + len(segment.syns) > self.maxSynapsesPerSegment) then
                !      numToFree = (len(segment.syns) + len(synsToAdd) - self.maxSynapsesPerSegment)
                !      segment.freeNSynapses(numToFree, inactiveSynIndices, self.verbosity)
                ! end if
                print '("WARNING: adapt_segment: create new synapses not implemented yet")'
            end if
        else ! Create a new segment
            block 
                integer(INT32) :: idx
                idx = columns(column) % dd_segment_count
                if (idx < N_DD_SEGMENTS_MAX) then
                    call inc(idx)
                    columns(column) % dd_segment_count = idx
                    columns(column) % sequence_segments(idx)        = seg_update % sequence_segment
                    columns(column) % dd_segment_destination(idx)   = seg_update % cell
                    columns(column) % dd_synapse_permanence(:, idx) = seg_update % dd_synapse_permanence
                    columns(column) % dd_synapse_origin(:, idx)     = seg_update % dd_synapse_origin
                else
                    print '("INFO:TP:adapt_segment: column"(1X,I0)": could not add segment, full")', column
                end if
            end block
        end if
        !if self.verbosity >= 4:
        !  print "   after:",
        !  segment.debugPrint()
    end subroutine adapt_segment
!========================================================================
end module htm_v3_tp
