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
module htm_v1_constants
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    use IFPORT, only: IRAND, RAND
 
    implicit none
    !=====================================================
        character(len=*),   parameter :: input_filename         = '../../Misc/data/AAAX_16x16/input.txt'
        integer(INT32),     parameter :: SENSOR_DIM1             = 16
        integer(INT32),     parameter :: SENSOR_DIM2             = 16
        integer(INT32),     parameter :: N_SENSORS              = SENSOR_DIM1 * SENSOR_DIM2
        integer(INT32),     parameter :: N_TIME_STEPS           = 200
        
        !=====================================================
        ! type constants
        integer(INT32),     parameter :: COLUMN_T               = INT32
        integer(INT32),     parameter :: CELL_ID_T              = INT32
        
        
        !=====================================================
        ! Spatial Pooler Parameters
        
        !defaultSpatialParams.put(KEY.INPUT_DIMENSIONS, new int[]{64});
        !defaultSpatialParams.put(KEY.POTENTIAL_RADIUS, 16);
        
        !defaultSpatialParams.put(KEY.INHIBITION_RADIUS, 0);
        !defaultSpatialParams.put(KEY.LOCAL_AREA_DENSITY, -1.0);
        
        !defaultSpatialParams.put(KEY.SYN_PERM_BELOW_STIMULUS_INC, 0.01);
        !defaultSpatialParams.put(KEY.SYN_PERM_TRIM_THRESHOLD, 0.5);
        !defaultSpatialParams.put(KEY.MIN_PCT_OVERLAP_DUTY_CYCLE, 0.001);
        !defaultSpatialParams.put(KEY.MIN_PCT_ACTIVE_DUTY_CYCLE, 0.001);

        ! The period used to calculate duty cycles. Higher values make it take longer to respond to changes in
        ! boost or synPerConnectedCell. Shorter values make it more unstable and likely to oscillate.
        !defaultSpatialParams.put(KEY.DUTY_CYCLE_PERIOD, 1000);
        
        
        
        
        
        ! If true, then during inhibition phase the winning columns are selected as the most active columns from
        ! the region as a whole. Otherwise, the winning columns are selected with respect to their local neighborhoods.
        logical(1),         parameter :: GLOBAL_INHIBITION      = .TRUE.
        logical(1),         parameter :: LEARNING_ON            = .TRUE.
        
        integer(INT32),     parameter :: COLUMN_DIM_1           = 100!2048
        integer(INT32),     parameter :: N_CELLS_PC     = 32     ! number of cell in a column
        integer(INT32),     parameter :: NUM_ACTIVE_COLUMNS_PER_INH_AREA = 40    ! number of columns that end up begin active within its inhibition radius

        ! This is a number specifying the minimum number of synapses that must be on in order for a columns to
        ! turn ON. The purpose of this is to prevent noise input from activating columns. Specified as a percent
        ! of a fully grown synapse.
        !defaultSpatialParams.put(KEY.STIMULUS_THRESHOLD, 0.0);
        
        ! The percent of the inputs, within a column's potential radius, that a column can be connected to.
        ! If set to 1, the column will be connected to every input within its potential radius. This parameter is
        ! used to give each column a unique potential pool when a large potentialRadius causes overlap between the
        ! columns. 
        real(REAL32),       parameter :: POTENTIAL_PCT          = 0.8

        ! Any synapse whose permanence value is above the connected threshold is a "connected synapse", meaning it can contribute to the cell's firing.
        real(REAL32),       parameter :: SYN_PERM_CONNECTED     = 0.1

        ! The amount by which an active synapse is incremented in each round. Specified as a percent of a fully grown synapse.
        real(REAL32),       parameter :: SYN_PERM_ACTIVE_INC    = 0.0001

        ! The amount by which an inactive synapse is decremented in each round. Specified as a percent of a fully grown synapse.
        real(REAL32),       parameter :: SYN_PERM_INACTIVE_DEC  = 0.0005

        ! The maximum overlap boost factor. Each column's overlap gets multiplied by a boost factor
        ! before it gets considered for inhibition. The actual boost factor for a column is number
        ! between 1.0 and maxBoost.
        real(REAL32),       parameter :: MAX_BOOST              = 1.0

        !=====================================================
        ! Temporal Memory Parameters
        !defaultTemporalParams.put(KEY.COLUMN_DIMENSIONS, new int[]{2048});
        !defaultTemporalParams.put(KEY.CELLS_PER_COLUMN, 32);
        
        !defaultTemporalParams.put(KEY.LEARNING_RADIUS, 2048);
        
        ! The maximum number of synapses added to a segment during learning.
        !defaultTemporalParams.put(KEY.MAX_NEW_SYNAPSE_COUNT, 20);
        
        ! Amount by which active permanences of synapses of previously predicted but inactive segments are decremented.
        !defaultTemporalParams.put(KEY.PREDICTED_SEGMENT_DECREMENT, 0.0);

        integer(INT32),     parameter :: MAX_NEW_SYNAPSE_COUNT  = 20
        
        ! Amount by which permanence of synapses are incremented during learning.
        real(REAL32),       parameter :: PERMANENCE_INCREMENT   = 0.01
        
        ! Amount by which permanences of synapses are decremented during learning.
        real(REAL32),       parameter :: PERMANENCE_DECREMENT   = 0.01
        
        ! If the number of synapses active on a segment is at least this threshold, it is selected as the best matching cell in a bursting column.
        integer(INT32),     parameter :: MIN_THRESHOLD          = 9

        ! If the number of active connected synapses on a segment is at least this threshold, the segment is said to be active.
        ! minimum overlap of connected synapses with active inputs 
        integer(INT32),     parameter :: PROXIMAL_DENDRITE_SEGMENT_ACTIVATION_THRESHOLD = 6
        integer(INT32),     parameter :: DISTAL_DENDRITE_SEGMENT_ACTIVATION_THRESHOLD   = 6
        
        ! Initial permanence of a new synapse
        real(REAL32),       parameter :: INITIAL_PERMANENCE         = 0.25

        ! if the permanence value for a synapse is greater than this value, it is said to be connected
        real(REAL32),       parameter :: CONNECTED_PERMANENCE       = 0.2 !NINT(0.2 * 256, INT8) - 127
    
        !=====================================================
        
        integer(INT32),     parameter :: N_DD_SYNAPSES          = 50    ! number of potential synapses per distal dendrite
        integer(INT32),     parameter :: N_PD_SYNAPSES          = 100    ! number of potential synapses per proximal dendrite
        
        integer(INT32),     parameter :: N_DD_SEGMENTS_MAX      = 1000   ! max number of distal dendrite segments available in a column
        integer(INT32),     parameter :: N_DD_SEGMENTS_TMP_MAX  = 100     ! max number of distal dendrite segments tmp available in a column
        
        integer(INT32),     parameter :: N_COLUMNS              = COLUMN_DIM_1    ! number columns in a region (layer)
        integer(INT32),     parameter :: N_CELLS                = N_CELLS_PC * N_COLUMNS    ! number of cell in a region (region is a layer)
        
    
        logical(1)          :: global_sensor_activity(N_SENSORS)        = .FALSE.
        integer(INT32)      :: global_inferred_sensor_activity(N_SENSORS)  = 0
        
        ! the following tho datastructures are send to each process
        logical(1)          :: global_column_is_active(N_COLUMNS)             = .FALSE.
        !integer(INT32)     :: global_cell_is_active_curr(N_COLUMNS)
        logical(1)          :: global_cell_is_active_curr(N_CELLS_PC, N_COLUMNS) = .FALSE.
        logical(1)          :: global_cell_is_active_prev(N_CELLS_PC, N_COLUMNS) = .FALSE.
        
        type distal_dendrite_segment_tmp_t
            integer(INT8)   :: cell                         = -1    ! cell to which this segment belongs, <0 if segment is invalid
            logical(1)      :: is_sequence                  = .FALSE.
            integer(INT32)  :: idx                          = -1    ! index of existing segment; -1 for new segments
            real(REAL32)    :: permanance(N_DD_SYNAPSES)    = -1
            integer(INT32)  :: origin(N_DD_SYNAPSES)        = -1
        end type distal_dendrite_segment_tmp_t
        
        type column_t
            !integer(INT32) :: predictive_state_curr, predictive_state_prev  = 0
            !integer(INT32) :: active_state_curr, active_state_prev          = 0
            !integer(INT32) :: learn_state_curr, learn_state_prev            = 0
            logical(1)      :: predictive_state_curr(N_CELLS_PC)    = .FALSE.
            logical(1)      :: predictive_state_prev(N_CELLS_PC)    = .FALSE.
            !logical(1)      :: active_state_curr(N_CELLS_PC)        = .FALSE.
            !logical(1)      :: active_state_prev(N_CELLS_PC)        = .FALSE.
            logical(1)      :: learn_state_curr(N_CELLS_PC)         = .FALSE.
            logical(1)      :: learn_state_prev(N_CELLS_PC)         = .FALSE.

            real(REAL32)    :: boost    = 1.0
            integer(INT32)  :: overlap  = 0

            real(REAL32)    :: proximal_dendrite_synapse_permanance( N_PD_SYNAPSES) = -128
            integer(INT32)  :: proximal_dendrite_synapse_origin(     N_PD_SYNAPSES) = -1

            integer(INT32)  :: distal_dendrite_segment_size                                 = 0
            logical(1)      :: distal_dendrite_segment_active_curr(   N_DD_SEGMENTS_MAX)    = .FALSE.
            logical(1)      :: distal_dendrite_segment_active_prev(   N_DD_SEGMENTS_MAX)    = .FALSE.
            integer(INT16)  :: distal_dendrite_segment_activity_curr( N_DD_SEGMENTS_MAX)    = 0 ! number of active synapses in this segment
            integer(INT16)  :: distal_dendrite_segment_activity_prev( N_DD_SEGMENTS_MAX)    = 0
            integer(INT8)   :: distal_dendrite_segment_destination(   N_DD_SEGMENTS_MAX)    = -1 ! cell to which this segment belongs to
           
            logical(1)      :: distal_dendrite_synapse_connected(     N_DD_SYNAPSES, N_DD_SEGMENTS_MAX) = .FALSE.
            real(REAL32)    :: distal_dendrite_synapse_permanance(    N_DD_SYNAPSES, N_DD_SEGMENTS_MAX) = -128
            integer(INT32)  :: distal_dendrite_synapse_origin(        N_DD_SYNAPSES, N_DD_SEGMENTS_MAX) = -1
            
            type(distal_dendrite_segment_tmp_t)  :: distal_dendrite_segment_tmp(N_DD_SEGMENTS_TMP_MAX)
            
            logical(1)      :: distal_dendrite_segment_sequence_curr(N_DD_SEGMENTS_MAX)                 = .FALSE.
            logical(1)      :: distal_dendrite_segment_sequence_prev(N_DD_SEGMENTS_MAX)                 = .FALSE.
            
        end type column_t
        type(column_t), allocatable :: columns(:)
        
        logical(1)      :: trace_on     = .TRUE.
        integer(INT32)  :: trace_cell   = 1
        integer(INT32)  :: trace_column = 1
        
    contains
!=====================================================
end module htm_v1_constants