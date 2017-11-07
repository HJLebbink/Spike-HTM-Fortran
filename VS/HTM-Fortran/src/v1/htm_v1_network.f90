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
module htm_v1_network
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use htm_v1_constants
    use htm_v1_tools
    use htm_v1_encoder, only: fetch_sensor_data, encode_pass_through
    use htm_v1_sp, only: spatial_pooling, update_inferred_sensor_activity
    use htm_v1_tp, only: temporal_pooling

    implicit none
    contains
    
!=====================================================
    subroutine load_input()
        call encode_pass_through(input_filename, SENSOR_DIM1, SENSOR_DIM2)
    end subroutine
!=====================================================
    subroutine run()
        integer(INT32) :: t, column_i
        
        do t = 1,N_TIME_STEPS
            global_cell_is_active_prev = global_cell_is_active_curr
            global_cell_is_active_curr = .FALSE.
            
            do column_i = 1,N_COLUMNS
                call advance_time(column_i)
            end do

            call fetch_sensor_data(t, global_sensor_activity)
            
            call spatial_pooling()
            call temporal_pooling()
            
            !call update_inferred_sensor_activity()
!            call to_string()
        end do
    end subroutine run
!=====================================================
    subroutine init()
        integer(INT32) :: column_i, cell_i, synapse_i
        integer(INT8)  :: allocateStatus
        
        allocate(columns(1:N_COLUMNS), STAT = allocateStatus)
        if (allocateStatus /= 0) stop "*** Not enough memory ***"
        
        ! init permanence values
        do column_i = 1,N_COLUMNS
            do synapse_i = 1,N_PD_SYNAPSES
                columns(column_i) % proximal_dendrite_synapse_permanance(synapse_i) = INITIAL_PERMANENCE
                columns(column_i) % proximal_dendrite_synapse_origin(synapse_i) = rand_int32(1, N_SENSORS)
            end do

            if (.FALSE.) then
                ! add one random distal dendrite segment
                columns(column_i) % distal_dendrite_segment_size = 1
            
                do cell_i = 1,N_CELLS_PC
                    columns(column_i) % distal_dendrite_segment_destination(1) = cell_i
                    do synapse_i = 1,N_DD_SYNAPSES
                        columns(column_i) % distal_dendrite_synapse_origin(    synapse_i, 1) = get_random_cell()
                        columns(column_i) % distal_dendrite_synapse_permanance(synapse_i, 1) = INITIAL_PERMANENCE
                        columns(column_i) % distal_dendrite_synapse_connected( synapse_i, 1) = (columns(column_i) % distal_dendrite_synapse_permanance(synapse_i, 1) > CONNECTED_PERMANENCE)
                    end do
                end do
            end if
        end do
    end subroutine
!=====================================================
    subroutine advance_time(column)
        integer(INT32), intent(in) :: column
        
        columns(column) % predictive_state_prev     = columns(column) % predictive_state_curr ! the only place were predictive state prev is written
        columns(column) % predictive_state_curr     = .FALSE.
        
        columns(column) % learn_state_prev          = columns(column) % learn_state_curr
        columns(column) % learn_state_curr          = .FALSE.
        
        columns(column) % distal_dendrite_segment_sequence_prev     = columns(column) % distal_dendrite_segment_sequence_curr
        columns(column) % distal_dendrite_segment_sequence_curr     = .FALSE.
                          
        columns(column) % distal_dendrite_segment_activity_prev     = columns(column) % distal_dendrite_segment_activity_curr
       !columns(column) % distal_dendrite_segment_activity_curr       is updated by update_distal_dendrite_segment_activity
        
        columns(column) % distal_dendrite_segment_active_prev       = columns(column) % distal_dendrite_segment_active_curr
       !columns(column) % distal_dendrite_segment_active_curr       = .FALSE. NOT needed
        
    end subroutine advance_time
!=====================================================
end module htm_v1_network