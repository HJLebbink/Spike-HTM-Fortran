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
module htm_v3_sp
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v3_constants
    use htm_v3_tools
    use htm_v3_print
    implicit none
    
        real(REAL32) :: active_duty_cycles          (N_COLUMNS) = 0.0
        real(REAL32) :: min_active_duty_cycles      (N_COLUMNS) = 0.0
        
        real(REAL32) :: overlap_duty_cycles         (N_COLUMNS) = 0.0
        real(REAL32) :: min_overlap_duty_cycles     (N_COLUMNS) = 0.0
        
        ! duty_cycle_period The period used to calculate duty cycles.
        ! Higher values make it take longer to respond to changes in
        ! boost. Shorter values make it potentially more unstable and
        ! likely to oscillate.
        integer(INT32), parameter :: DUTY_CYCLE_PERIOD      = 1000

        real(REAL32),   parameter :: MAX_BOOST              = 5.0
        
        ! minimum tolerated activity duty cycle, given as percent of neighbors' 
        ! activity duty cycle. minPctActiveDutyCycles real number of the 
        ! minimum tolerated activity duty cycle.
        real(REAL32),   parameter :: MIN_PCT_ACTIVE_DUTY_CYCLES = 0.1!0.001
        
        ! minimum tolerated overlaps, given as percent of neighbors overlap 
        ! score. minPctOverlapDutyCycles real number of the minimum tolerated 
        ! overlaps.
        real(REAL32),   parameter :: MIN_PCT_OVERLAP_DUTY_CYCLES = 0.15!0.001

        integer(INT32), parameter :: UPDATE_PERIOD          = 1!50
        
        integer(INT32) :: iteration_num                     = 0
        integer(INT32) :: iteration_learn_num               = 0
        
        ! Permanence increment amount for columns that have not been recently active.
        ! real number of the permanence increment amount for columns that have not been
        ! recently active, must be larger than 0.
        real(REAL32), parameter :: SYN_PERM_BELOW_STIMULS_INC = 0.01
        
    contains
!========================================================================
    subroutine compute_sp(                                              &
            sensor_activity,                                            &
            learn,                                                      &
            ! out
            active_columns)
            
    !This is the main workshorse method of the SpatialPooler class. This
    !method takes an input vector and computes the set of output active
    !columns. If 'learn' is set to True, this method also performs
    !learning.

    !@param inputVector An array of integer 0's and 1's that comprises
    !    the input to the spatial pooler. The length of the
    !    array must match the total number of input bits implied by
    !    the constructor (also returned by the method getNumInputs). In
    !    cases where the input is multi-dimensional, inputVector is a
    !    flattened array of inputs.

    !@param learn A boolean value indicating whether learning should be
    !    performed. Learning entails updating the permanence values of
    !    the synapses, duty cycles, etc. Learning is typically on but
    !    setting learning to 'off' is useful for analyzing the current
    !    state of the SP. For example, you might want to feed in various
    !    inputs and examine the resulting SDR's. Note that if learning
    !    is off, boosting is turned off and columns that have never won
    !    will be removed from activeVector.  TODO: we may want to keep
    !    boosting on even when learning is off.

    !@param activeVector An array representing the winning columns after
    !    inhibition. The size of the array is equal to the number of
    !    columns (also returned by the method getNumColumns). This array
    !    will be populated with 1's at the indices of the active columns,
    !    and 0's everywhere else. In the case where the output is
    !    multi-dimensional, activeVector represents a flattened array
    !    of outputs.

        logical(1), intent(in)  :: sensor_activity  (N_SENSORS)
        logical(1), intent(in)  :: learn
        logical(1), intent(out) :: active_columns   (N_COLUMNS)
        
        ! local variables
        integer(INT32) :: column_i
        integer(INT32) :: overlap                   (N_COLUMNS)
        real(REAL32)   :: boosted_overlap           (N_COLUMNS)
        
        iteration_num = iteration_num + 1
        if (learn) iteration_learn_num = iteration_learn_num + 1
        
        do column_i = 1, N_COLUMNS
            call calc_overlap(                                          &
                column              = column_i,                         &
                sensor_activity     = sensor_activity,                  &
                overlap             = overlap(column_i))
            
            if (learn) then 
                boosted_overlap(column_i) = overlap(column_i) + (rand() * 0.001)
                !boosted_overlap(column_i) = (overlap(column_i) * columns(column_i) % boost_factor) + (rand() * 0.001) ! add a small rand number to allow the selection of the top n best overlaps to be done in parallel
            else
                boosted_overlap(column_i) = overlap(column_i) + (rand() * 0.001)
            end if
        end do
        
        ! barrier: the overlap of all columns needs to be known
        
        do column_i = 1, N_COLUMNS
            call inhibition(                                            &
                column              = column_i,                         &
                learn               = learn,                            &
                boosted_overlap     = boosted_overlap,                  &
                sensor_activity     = sensor_activity,                  &
                ! out
                active_column       = active_columns(column_i))

            if (learn) then
                call update_duty_cycles(column_i, overlap(column_i), active_columns(column_i))
                call bump_up_weak_columns(column_i, overlap_duty_cycles(column_i), min_overlap_duty_cycles(column_i))
                call update_boost_factors(column_i, active_duty_cycles(column_i), min_active_duty_cycles(column_i))
                if (is_update_round()) then
                    call update_inhibition_radius(column_i)
                    call update_min_duty_cycles(column_i)
                end if
            end if
        end do
#       if _DEBUG
            if (.TRUE.) then
                write(*, '("INFO:SP:compute_sp: sensor activity IN:")', advance='yes')
                call print_sensor_activity(sensor_activity, SENSOR_DIM1, SENSOR_DIM2)
                write(*, '("INFO:SP:compute_sp: active columns OUT:")', advance='no')
                call print_bitset(active_columns)
            end if
#       endif
    end subroutine compute_sp
!========================================================================
    subroutine calc_overlap(                                            &
            column,                                                     &
            sensor_activity,                                            &
            ! out
            overlap)
            
        integer(INT32), intent(in)  :: column
        logical(1),     intent(in)  :: sensor_activity(N_SENSORS)
        integer(INT32), intent(out) :: overlap
        ! local variables
        integer(INT32) :: synapse_i, sensor_idx
        
        if (.FALSE.) then
            overlap = 0
            do synapse_i = 1, N_PD_SYNAPSES
                if (is_PD_connected(column, synapse_i, SP_PD_PERMANENCE_THRESHOLD)) then
                    sensor_idx = columns(column) % pd_synapse_origin(synapse_i)
                    if (sensor_activity(sensor_idx)) then
                        overlap = overlap + 1
                    end if
                end if
            end do
        else
            block
                integer(INT32) :: x(N_PD_SYNAPSES)
                integer(INT32) :: overlap2
                
                x = (columns(column) % pd_synapse_permanence >= SP_PD_PERMANENCE_THRESHOLD).AND.sensor_activity(columns(column) % pd_synapse_origin)
                overlap = SUM(x)
                
#               if _DEBUG
                overlap2 = 0
                do synapse_i = 1, N_PD_SYNAPSES
                    if (is_PD_connected(column, synapse_i, SP_PD_PERMANENCE_THRESHOLD)) then
                        sensor_idx = columns(column) % pd_synapse_origin(synapse_i)
                        if (sensor_activity(sensor_idx)) then
                            overlap2 = overlap2 + 1
                        end if
                    end if
                end do
                
                if (overlap /= overlap2) then
                    print '("ERROR:SP:calc_overlap: ERROR:")'
                    read *
                end if
#               endif

            end block
        end if
        
        if (overlap < STIMULUS_THRESHOLD) then
            overlap = 0
        end if
        
#       if _DEBUG
            if (.FALSE.) print '("INFO:SP:calc_overlap: column",(1X,I0)," has overlap = ",(1X,I0))', column, overlap
#       endif
    end subroutine calc_overlap
!========================================================================
    subroutine inhibition(                                              &
            column,                                                     &
            learn,                                                      &
            boosted_overlap,                                            &
            sensor_activity,                                            &
            ! out
            active_column)
            
        integer(INT32), intent(in)  :: column
        logical(1),     intent(in)  :: learn
        real(REAL32),   intent(in)  :: boosted_overlap      (N_COLUMNS)
        logical(1),     intent(in)  :: sensor_activity      (N_SENSORS)
        logical(1),     intent(out) :: active_column
        
        ! local variables
        integer(INT32) :: synapse_i, sensor_idx
        real(REAL32)   :: old_permanence, new_permanence
        logical(1)     :: is_inhibited
        
        active_column = .FALSE.

        if (boosted_overlap(column) > 0.0) then
            is_inhibited = column_is_inhibited(column, boosted_overlap, INHIBITION_TOP)
            if (is_inhibited) then
                !print '("INFO:SP:inhibition: column ",i4," is inhibited")', column
            else 
#               if _DEBUG
                if (.FALSE.) then
                    print '("INFO:SP:inhibition: column",(1X,I0)," is not inhibited; has most (overlap=",f,"; INHIBITION_TOP=",(1X,I0),").")', column, boosted_overlap(column), INHIBITION_TOP
                end if
#               endif
                
                active_column = .TRUE.

                if (learn) then
                    do synapse_i = 1,N_PD_SYNAPSES
                        sensor_idx = columns(column) % pd_synapse_origin(synapse_i)
                        old_permanence = columns(column) % pd_synapse_permanence(synapse_i)

                        !print '("INFO:SP:inhibition:      old permanance for synpase ",i5," in column ",i4,":",f,").")', synapse_i, column_i, columns(column_i) % pd_synapse_permanance(synapse_i)
                        if (sensor_activity(sensor_idx)) then
                            new_permanence = MIN(1.0, old_permanence + PD_PERMANENCE_INC)
                            !print '("INFO:SP:inhibition: inc: new permanance for synpase ",i5," in column ",i4,":",f,").")', synapse_i, column, new_permanence
                        else
                            new_permanence = MAX(0.0, old_permanence - PD_PERMANENCE_DEC)
                            !print '("INFO:SP:inhibition: dec: new permanance for synpase ",i5," in column ",i4,":",f,").")', synapse_i, column, new_permanence
                        end if
                                
#                       if _DEBUG
                            if (.FALSE.) then
                                if (old_permanence > SP_PD_PERMANENCE_THRESHOLD) then
                                    if (new_permanence <= SP_PD_PERMANENCE_THRESHOLD) then
                                        print '("INFO:SP:inhibition: column ",i3, ", synapse ",i4," went from connected to not connected")', column, synapse_i
                                    end if
                                else
                                    if (new_permanence > SP_PD_PERMANENCE_THRESHOLD) then
                                        print '("INFO:SP:inhibition: column ",i3, ", synapse ",i4," went from not connected to connected")', column, synapse_i
                                    end if
                                end if
                            end if
#                       endif
                        columns(column) % pd_synapse_permanence(synapse_i) = new_permanence
                    end do
                end if
            end if
        end if
            
        ! TODO: 1] update network_boost(column_i)
        ! TODO: 2] increase permanence if input overlap is small
    end subroutine inhibition
!========================================================================
    function column_is_inhibited(                                       &
            column,                                                     &
            boosted_overlap,                                            &
            inhibition_top) result(r)
            
        integer(INT32), intent(in) :: column
        real(REAL32),   intent(in) :: boosted_overlap       (N_COLUMNS)
        integer(INT32), intent(in) :: inhibition_top
        logical(1) :: r
        
        integer(INT32) :: n_columns_with_more_overlap
        real(REAL32)   :: overlap_this_column
        integer(INT32) :: tmp_sum                   (N_COLUMNS)
        
        overlap_this_column = boosted_overlap(column)
        
        if (GLOBAL_INHIBITION) then ! this can be done faster by just returing the k highest overlaps
            tmp_sum = (boosted_overlap(:) > overlap_this_column)
            n_columns_with_more_overlap = SUM(tmp_sum)
            r = (n_columns_with_more_overlap >= inhibition_top) 
            !print '("INFO: column_is_inhibited: column ",i4," is inhibited ",i1, "; n_columns_with_more_overlap=",i4)', column, r, n_columns_with_more_overlap
        else
            !1] retrieve the columns in the neigborhood of column_i
            r = .TRUE.
            print '("not implemented yet")'
            read *
        end if
    end function column_is_inhibited
!========================================================================
    elemental function is_PD_connected(column, synapse, permanence_threshold) result(r)
        integer(INT32), intent(in) :: column, synapse
        real(REAL32),   intent(in) :: permanence_threshold
        logical(1) :: r
        r = (columns(column) % pd_synapse_permanence(synapse) >= permanence_threshold)
    end function is_PD_connected
!========================================================================
    subroutine update_duty_cycles(                                      &
            column,                                                     &
            overlap,                                                    &
            active_column)
    
    ! Updates the duty cycles for each column. The OVERLAP duty cycle is a moving
    ! average of the number of inputs which overlapped with the each column. The
    ! ACTIVITY duty cycles is a moving average of the frequency of activation for
    ! each column.

    ! overlap. an int vector containing the overlap score for each column.
    ! The overlap score for a column is defined as the number
    ! of synapses in a "connected state" (connected synapses)
    ! that are connected to input bits which are turned on.

    ! active_columns.  An bit array containing the the active columns,
    ! the set of columns which survived inhibition

        integer(INT32), intent(in) :: column
        integer(INT32), intent(in) :: overlap
        logical(1),     intent(in) :: active_column

    ! vector<UInt> newOverlapVal(numColumns_, 0);
    ! vector<UInt> newActiveVal(numColumns_, 0);

    ! for (UInt i = 0; i < numColumns_; i++) {
    !   newOverlapVal[i] = overlaps[i] > 0 ? 1 : 0;
    !   newActiveVal[i] = activeArray[i] > 0 ? 1 : 0;
    ! }
    ! UInt period = dutyCyclePeriod_ > iterationNum_ ? iterationNum_ : dutyCyclePeriod_;
    ! updateDutyCyclesHelper_(overlapDutyCycles_, newOverlapVal, period);
    ! updateDutyCyclesHelper_(activeDutyCycles_, newActiveVal, period);
    !void SpatialPooler::updateDutyCyclesHelper_(vector<Real>& dutyCycles, vector<UInt>& newValues, UInt period) {
    !   for (UInt i = 0; i < dutyCycles.size(); i++) {
    !       dutyCycles[i] = (dutyCycles[i] * (period - 1) + newValues[i]) / period;
    !   }
    !}
        integer(INT32) :: period
        period = MIN(iteration_num, DUTY_CYCLE_PERIOD)

        if (overlap > 0) then
            overlap_duty_cycles(column) = (overlap_duty_cycles(column) * (period - 1) + 1) / period
        else
            overlap_duty_cycles(column) = (overlap_duty_cycles(column) * (period - 1) + 0) / period
        end if
            
        if (active_column) then
            active_duty_cycles(column)  = (active_duty_cycles(column)  * (period - 1) + 1) / period
        else
            active_duty_cycles(column)  = (active_duty_cycles(column)  * (period - 1) + 0) / period
        end if
    end subroutine update_duty_cycles
!========================================================================
    subroutine update_boost_factors(column, active_duty_cycle, min_active_duty_cycle)
        integer(INT32), intent(in) :: column
        real(REAL32),   intent(in) :: active_duty_cycle
        real(REAL32),   intent(in) :: min_active_duty_cycle

    ! Update the boost factors for all columns. The boost factors are used to
    ! increase the overlap of inactive columns to improve their chances of
    ! becoming active. and hence encourage participation of more columns in the
    ! learning process. This is a line defined as: y = mx + b boost =
    ! (1-maxBoost)/minDuty * dutyCycle + maxFiringBoost. Intuitively this means
    ! that columns that have been active enough have a boost factor of 1, meaning
    ! their overlap is not boosted. Columns whose active duty cycle drops too much
    ! below that of their neighbors are boosted depending on how infrequently they
    ! have been active. The more infrequent, the more they are boosted. The exact
    ! boost factor is linearly interpolated between the points (dutyCycle:0,
    ! boost:maxFiringBoost) and (dutyCycle:minDuty, boost:1.0).
    !
    !         boostFactor
    !             ^
    ! maxBoost _  |
    !             |\
    !             | \
    !       1  _  |  \ _ _ _ _ _ _ _
    !             |
    !             +--------------------> activeDutyCycle
    !                |
    !         minActiveDutyCycle
   
!    for (UInt i = 0; i < numColumns_; i++) {
!        if (minActiveDutyCycles_[i] <= 0) {
!            continue;
!        }
!        if (activeDutyCycles_[i] > minActiveDutyCycles_[i]) {
!            boostFactors_[i] = 1.0;
!            continue;
!        }
!        boostFactors_[i] = ((1 - maxBoost_) / minActiveDutyCycles_[i] * activeDutyCycles_[i]) + maxBoost_;
!    }
        if (min_active_duty_cycle > 0) then

#           if _DEBUG
                if (.FALSE.) print '("INFO:SP:update_boost_factors: column",(1X,I0)," old boost factor ",f)', column, columns(column) % boost_factor
#           endif

            if (active_duty_cycle > min_active_duty_cycle) then
                columns(column) % boost_factor = 1.0
            else 
                columns(column) % boost_factor = ((1.0 - MAX_BOOST) / min_active_duty_cycle * active_duty_cycle) + MAX_BOOST
#               if _DEBUG
                    if (.FALSE.) print '("INFO:SP:update_boost_factors: column",(1X,I0)," new boost factor ",f,"; madc=",f,"; adc=",f)', column, columns(column) % boost_factor, min_active_duty_cycle, active_duty_cycle
#               endif
            end if
        end if
   end subroutine update_boost_factors
!========================================================================
    subroutine bump_up_weak_columns(column, overlap_duty_cycle, min_overlap_duty_cycle)
        integer(INT32), intent(in) :: column
        real(REAL32),   intent(in) :: overlap_duty_cycle
        real(REAL32),   intent(in) :: min_overlap_duty_cycle
        
     ! This method increases the permanence values of synapses of columns whose
     ! activity level has been too low. Such columns are identified by having an
     ! overlap duty cycle that drops too much below those of their peers. The
     ! permanence values for such columns are increased.
    
        integer(INT32) :: synapse_i
#       if _DEBUG
            if (.FALSE.) print '("INFO:SP:bump_up_weak_columns: column ="(1X,I0)"; overlap_duty_cycle=",f,"; min_overlap_duty_cycle=",f)', column, overlap_duty_cycle, min_overlap_duty_cycle
#       endif
        
        if (overlap_duty_cycle < min_overlap_duty_cycle) then
            
#           if _DEBUG
                if (.TRUE.) print '("INFO:SP:bump_up_weak_columns: column ="(1X,I0)" is bumped up; overlap_duty_cycle=",f,"; min_overlap_duty_cycle=",f)', column, overlap_duty_cycle, min_overlap_duty_cycle
#           endif
            
            do synapse_i = 1, N_PD_SYNAPSES
                if (NOT(is_PD_connected(column, synapse_i, SP_PD_PERMANENCE_THRESHOLD))) then
                    if (columns(column) % pd_synapse_origin(synapse_i) /= -1) then
                        columns(column) % pd_synapse_permanence(synapse_i) = &
                            columns(column) % pd_synapse_permanence(synapse_i) + PD_PERMANENCE_INC_WEAK
                    end if
                end if
            end do
        end if
    end subroutine bump_up_weak_columns
!========================================================================
    subroutine update_inhibition_radius(column)
        integer(INT32), intent(in) :: column
        !TODO
    end subroutine update_inhibition_radius
!========================================================================
    subroutine update_min_duty_cycles(column)
        integer(INT32), intent(in) :: column

        if (GLOBAL_INHIBITION) then !.OR.(inhibitionRadius > N_COLUMNS) then
            !Real maxActiveDutyCycles = *max_element(activeDutyCycles_.begin(), activeDutyCycles_.end());
            !Real maxOverlapDutyCycles = *max_element(overlapDutyCycles_.begin(), overlapDutyCycles_.end());
            !fill(minActiveDutyCycles_.begin(), minActiveDutyCycles_.end(), minPctActiveDutyCycles_ * maxActiveDutyCycles);
            !fill(minOverlapDutyCycles_.begin(), minOverlapDutyCycles_.end(), minPctOverlapDutyCycles_ * maxOverlapDutyCycles);
            
            min_active_duty_cycles(column)  = MIN_PCT_ACTIVE_DUTY_CYCLES  * MAXVAL(active_duty_cycles)
            min_overlap_duty_cycles(column) = MIN_PCT_OVERLAP_DUTY_CYCLES * MAXVAL(overlap_duty_cycles)
#           if _DEBUG
                if (.FALSE.) print '("INFO:SP:update_min_duty_cycles: min_active_duty_cycles = ",f,"; min_overlap_duty_cycles = ",f)', min_active_duty_cycles(1), min_overlap_duty_cycles(1)
#           endif
        else
            print '("not implemented yet")'
            read *
            
            !for (UInt i = 0; i < numColumns_; i++) {
            !vector<UInt> neighbors;

            !getNeighborsND_(i, columnDimensions_, inhibitionRadius_, false, neighbors);
            !neighbors.push_back(i);
            !Real maxActiveDuty = 0;
            !Real maxOverlapDuty = 0;
            !for (auto & neighbor : neighbors) {
                !UInt index = neighbor;
                !maxActiveDuty = max(maxActiveDuty, activeDutyCycles_[index]);
                !maxOverlapDuty = max(maxOverlapDuty, overlapDutyCycles_[index]);
            !}
            !minActiveDutyCycles_[i] = maxActiveDuty * minPctActiveDutyCycles_;
            !minOverlapDutyCycles_[i] = maxOverlapDuty * minPctOverlapDutyCycles_;

        end if
    end subroutine update_min_duty_cycles
!========================================================================
    function is_update_round() result(r)
        logical(1) :: r 
        r = (MODULO(iteration_num, UPDATE_PERIOD) == 0)
    end function is_update_round
!========================================================================
end module htm_v3_sp
