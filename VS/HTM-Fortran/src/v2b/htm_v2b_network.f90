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
module htm_v2b_network
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use tools
    use htm_v2b_constants
    use htm_v2b_tools
    use htm_v2b_encoder, only: encode_pass_through, fetch_sensor_activity, print_sensor_activity2, print_inferred_sensor_activity, encoder_cleanup
    use htm_v2b_sp,      only: compute_sp
    use htm_v2b_tm,      only: compute_tm
    use htm_v2b_print
    
    implicit none

        integer(INT64), allocatable      :: active_cells_all             (:)
        integer(INT64), allocatable      :: prev_active_cells_all        (:)
        integer(INT64), allocatable      :: winner_cells_all             (:)
        integer(INT64), allocatable      :: prev_winner_cells_all        (:)
        integer(INT64), allocatable      :: predictive_cells_all         (:)
        integer(INT64), allocatable      :: prev_predictive_cells_all    (:)
        
    contains
!========================================================================
    subroutine load_input()
        call encode_pass_through(input_filename, SENSOR_DIM1, SENSOR_DIM2)
    end subroutine
!========================================================================
    subroutine run(n_time_steps, progress)
        integer(INT32), intent(in) :: n_time_steps
        logical(1),     intent(in) :: progress
        ! local variables
        logical(1),     parameter :: learn = .TRUE.

        integer(INT32)  :: t
        logical(1)      :: sensor_activity              (N_SENSORS)
        logical(1)      :: active_columns               (N_COLUMNS)
        integer(INT32)  :: total_mismatch

        call clear_all_64_array(active_cells_all)
        call clear_all_64_array(winner_cells_all)
        call clear_all_64_array(predictive_cells_all)
        
        total_mismatch = 0
        
        do t = 1,n_time_steps
            
            prev_active_cells_all           = active_cells_all
            prev_winner_cells_all           = winner_cells_all
            prev_predictive_cells_all       = predictive_cells_all
            
            call fetch_sensor_activity(t, sensor_activity)
            
            call compute_sp(                                                &
                sensor_activity             = sensor_activity,              &
                learn                       = learn,                        &
                ! out
                active_columns              = active_columns)
            
            call compute_tm(                                                &
                active_columns              = active_columns,               &
                prev_active_cells_all       = prev_active_cells_all,        &
                prev_winner_cells_all       = prev_winner_cells_all,        &
                prev_predictive_cells_all   = prev_predictive_cells_all,    &
                learn                       = learn,                        &
                ! out
                active_cells_all            = active_cells_all,             &
                winner_cells_all            = winner_cells_all,             &
                predictive_cells_all        = predictive_cells_all)
            
#if _DEBUG
                if (.FALSE.) then
                    write(*, '("INFO:run: active columns at t ="(1X,I0)":")', advance='yes'), t
                    call print_active_columns(active_columns, NINT(SQRT(REAL(N_COLUMNS))), NINT(SQRT(REAL(N_COLUMNS))))
                end if
                if (.FALSE.) then
                    write(*, '("INFO:run: active cells at t ="(1X,I0)": ")', advance='yes'), t
                    call print_active_cells(active_cells_all)
                end if
                if (.FALSE.) then
                    write(*, '("INFO:run: predictive cells at t ="(1X,I0)": ")', advance='yes'), t
                    call print_active_cells(predictive_cells_all)
                end if
#endif

            if (.TRUE.) then
                total_mismatch = total_mismatch + calc_mismatch(t, predictive_cells_all)
                if (MOD(t, 100) == 0) then
                    write(*, '("INFO:run: total mismatch at t ="(1X,I0)":"(1X,I0))', advance='yes'), t, total_mismatch
                    total_mismatch = 0
                end if
            end if

            if (progress) then
                call print_progress(t, active_cells_all, active_columns, predictive_cells_all)
            end if
        end do
        
        ! housekeeping
        call encoder_cleanup()
        call network_cleanup()
        
    end subroutine run
!========================================================================
    function calc_mismatch(t, predictive_cells) result(r)
        integer(INT32), intent(in) :: t
        integer(INT64), intent(in) :: predictive_cells      (N_COLUMNS)
        integer(INT32) :: r
        ! local variables
        integer(INT32)  :: sensor_activity_int              (N_SENSORS)
        logical(1)      :: sensor_activity1                 (N_SENSORS)
        logical(1)      :: sensor_activity2                 (N_SENSORS)
        
        call get_predicted_sensor_activity(predictive_cells, sensor_activity_int)
        sensor_activity1 = (sensor_activity_int > 0)
        call fetch_sensor_activity(t+1, sensor_activity2)
    
        sensor_activity_int = (sensor_activity1 /= sensor_activity2)
        r = SUM(sensor_activity_int)
    end function calc_mismatch
!========================================================================
    subroutine print_progress(t, active_cells, active_columns, predictive_cells)
        integer(INT32), intent(in) :: t
        integer(INT64), intent(in) :: active_cells          (N_COLUMNS)
        logical(1),     intent(in) :: active_columns        (N_COLUMNS)
        integer(INT64), intent(in) :: predictive_cells      (N_COLUMNS)
        ! local variables
        integer(INT32), parameter :: progress_frequency = 1
        logical(1),     parameter :: show_mismatch      = .TRUE.
        
        integer(INT32)  :: inferred_sensor_activity     (N_SENSORS)
        logical(1)      :: sensor_activity1             (N_SENSORS)
        logical(1)      :: sensor_activity2             (N_SENSORS)
        
        !if (t > 6990) then
        if (MODULO(t, progress_frequency) == 0) then
            print '("=============================================================")'

            if (.TRUE.) then ! print active cells
                !print '("=====")'
                write(*, '("INFO:at t ="(1X,I0)": active cells:")', advance='yes'), t
                call print_active_cells(active_cells)
            end if
            
            if (.FALSE.) then ! print infered sensor activity
                print '("=====")'
                write(*, '("INFO:at t ="(1X,I0)": inferred sensor activity:")', advance='yes'), t
                call get_inferred_sensor_activity(active_columns, inferred_sensor_activity)
                sensor_activity1 = inferred_sensor_activity>0
                call print_sensor_activity(sensor_activity1, SENSOR_DIM1, SENSOR_DIM2)
            end if

            if (.TRUE.) then ! print predicted sensor activity
                !print '("=====")'
                print '("INFO:at t ="(1X,I0)": predicted sensor activity at t ="(1X,I0)":")', t, t+1
                print '("INFO:predicted activity | correct activity | mismatch")'
                call get_predicted_sensor_activity(predictive_cells, inferred_sensor_activity)
                sensor_activity1 = inferred_sensor_activity>0
                call fetch_sensor_activity(t+1, sensor_activity2)
                call print_sensor_activity2(sensor_activity1, sensor_activity2, SENSOR_DIM1, SENSOR_DIM2)
            end if
            
            if (.FALSE.) then
                block
                    real(REAL32)   :: boost_factors               (N_COLUMNS)
                    real(REAL32)   :: projected_boost_factors     (N_SENSORS)
                    integer(INT32) :: column_i
                    
                    do column_i = 1,N_COLUMNS
                        boost_factors(column_i) = columns(column_i) % boost_factor
                    end do
                    
                    call print_boost_factors(boost_factors)
                    call get_projected_boost_factors(boost_factors, projected_boost_factors)
                    call print_projected_boost_factors(projected_boost_factors, SENSOR_DIM1, SENSOR_DIM2)
                end block
            end if
        end if
    end subroutine print_progress
!========================================================================
    subroutine init(param)
        type(param_t), intent(in) :: param
        ! local variables
        integer(INT32) :: column_i, synapse_i
        integer(INT8)  :: allocateStatus
        
        allocate(columns(param % n_columns), STAT = allocateStatus)
        if (allocateStatus /= 0) stop "*** Not enough memory ***"

        allocate(active_cells_all(param % n_columns), STAT = allocateStatus)
        if (allocateStatus /= 0) stop "*** Not enough memory ***"
        allocate(prev_active_cells_all(param % n_columns), STAT = allocateStatus)
        if (allocateStatus /= 0) stop "*** Not enough memory ***"
        allocate(winner_cells_all(param % n_columns), STAT = allocateStatus)
        if (allocateStatus /= 0) stop "*** Not enough memory ***"
        allocate(prev_winner_cells_all(param % n_columns), STAT = allocateStatus)
        if (allocateStatus /= 0) stop "*** Not enough memory ***"
        allocate(predictive_cells_all(param % n_columns), STAT = allocateStatus)
        if (allocateStatus /= 0) stop "*** Not enough memory ***"
        allocate(prev_predictive_cells_all(param % n_columns), STAT = allocateStatus)
        if (allocateStatus /= 0) stop "*** Not enough memory ***"

        ! init permanence values
        do column_i = 1,param % n_columns
            columns(column_i) % distal_dendrite_segment_size = 0

            do synapse_i = 1,N_PD_SYNAPSES
                columns(column_i) % proximal_dendrite_synapse_permanence(synapse_i) = SP_PD_PERMANENCE_INIT
                columns(column_i) % proximal_dendrite_synapse_origin(synapse_i) = rand_int32(1, N_SENSORS)
            end do
            
            if (.FALSE.) then ! sort the proximal dendrites
                call sort(columns(column_i) % proximal_dendrite_synapse_origin)
            end if
        end do
    end subroutine init
!========================================================================
    subroutine get_inferred_sensor_activity(active_columns, inferred_sensor_activity)
        logical(1),     intent(in)  :: active_columns               (N_COLUMNS)
        integer(INT32), intent(out) :: inferred_sensor_activity     (N_SENSORS)
        ! local variables
        integer(INT32) :: column_i, synapse_i, sensor_idx
        
        inferred_sensor_activity = 0
        
        do column_i = 1,N_COLUMNS
            if (active_columns(column_i)) then
                do synapse_i = 1, N_PD_SYNAPSES
                    if (columns(column_i) % proximal_dendrite_synapse_permanence(synapse_i) >= SP_PD_PERMANENCE_THRESHOLD) then
                        sensor_idx = columns(column_i) % proximal_dendrite_synapse_origin(synapse_i)
                        inferred_sensor_activity(sensor_idx) = inferred_sensor_activity(sensor_idx) + 1
                    end if
                end do
            end if
        end do
    end subroutine get_inferred_sensor_activity
!========================================================================
    subroutine get_predicted_sensor_activity(predicted_cells_all, predicted_sensor_activity)
        integer(INT64), intent(in)  :: predicted_cells_all          (N_COLUMNS)
        integer(INT32), intent(out) :: predicted_sensor_activity    (N_SENSORS)
        ! local variables
        integer(INT32) :: column_i, synapse_i, sensor_idx
        integer(INT8)  :: cell_i
        integer(INT64) :: predicted_cells
        predicted_sensor_activity = 0
        
        do column_i = 1, N_COLUMNS
            predicted_cells = predicted_cells_all(column_i)
            if (is_not_empty_64(predicted_cells)) then 
                do cell_i = 1, N_CELLS_PC
                    if (get_element_64(predicted_cells, cell_i)) then
                        do synapse_i = 1, N_PD_SYNAPSES
                            if (columns(column_i) % proximal_dendrite_synapse_permanence(synapse_i) >= SP_PD_PERMANENCE_THRESHOLD) then
                                sensor_idx = columns(column_i) % proximal_dendrite_synapse_origin(synapse_i)
                                predicted_sensor_activity(sensor_idx) = predicted_sensor_activity(sensor_idx) + 1
                            end if
                        end do
                    end if
                end do
            end if
        end do
    end subroutine get_predicted_sensor_activity
!========================================================================
    subroutine get_projected_boost_factors(boost_factor, projected_boost_factor)
        real(REAL32), intent(in)  :: boost_factor           (N_COLUMNS)
        real(REAL32), intent(out) :: projected_boost_factor (N_SENSORS)
        ! local variable
        integer(INT32) :: synapse_i, column_i, sensor_idx, sensor_i
        integer(INT32) :: counter                           (N_SENSORS)
        projected_boost_factor = 0.001
        counter = 0
        
        do column_i = 1, N_COLUMNS
            do synapse_i = 1, N_PD_SYNAPSES
                if (columns(column_i) % proximal_dendrite_synapse_permanence(synapse_i) >= SP_PD_PERMANENCE_THRESHOLD) then
                    sensor_idx = columns(column_i) % proximal_dendrite_synapse_origin(synapse_i)
                    projected_boost_factor(sensor_idx) = projected_boost_factor(sensor_idx) + boost_factor(column_i)
                    counter(sensor_idx) = counter(sensor_idx) + 1
                end if
            end do
        end do
        
        do sensor_i = 1, N_SENSORS
            if (counter(sensor_i) > 0) then
                projected_boost_factor(sensor_i) = projected_boost_factor(sensor_i)/counter(sensor_i)
            end if
        end do
    end subroutine get_projected_boost_factors
!========================================================================
    subroutine network_cleanup()
        if (ALLOCATED(columns)) deallocate(columns)
        
        if (ALLOCATED(active_cells_all)) deallocate(active_cells_all)
        if (ALLOCATED(prev_active_cells_all)) deallocate(prev_active_cells_all)
        if (ALLOCATED(winner_cells_all)) deallocate(winner_cells_all)
        if (ALLOCATED(prev_winner_cells_all)) deallocate(prev_winner_cells_all)
        if (ALLOCATED(predictive_cells_all)) deallocate(predictive_cells_all)
        if (ALLOCATED(prev_predictive_cells_all)) deallocate(prev_predictive_cells_all)

    end subroutine network_cleanup
!========================================================================
    function estimate_memory_usage_byte(param) result(r)
        type(param_t), intent(in) :: param
        
        integer(INT64) :: r
        integer(INT64) :: column_size
        column_size =                                           &
            (3 * (1 * param % n_cells_pc)) +                    &
            (5 * (1 * N_DD_SEGMENTS_MAX)) +                     &
            (1 * 1) +                                           &
            ! proximal dendrites
            (2 * (4 * N_PD_SYNAPSES)) +                         &
            ! distal dendrites
            (1 * 4) +                                           &
            (1 * (4 * N_DD_SEGMENTS_MAX)) +                     &
            (1 * (4 * N_DD_SYNAPSES * N_DD_SEGMENTS_MAX)) +     &
            (1 * (4 * N_DD_SYNAPSES * N_DD_SEGMENTS_MAX))
        
        !print '("INFO: estimated column size"(1X,I0)" kbyte")', column_size/1024
        r = param % n_columns * column_size 
    end function 
!========================================================================
end module htm_v2b_network
