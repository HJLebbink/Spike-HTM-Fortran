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
module htm_v3_constants
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT8, INT16, INT32, INT64, REAL32, REAL64
    implicit none

    !========================================================================
    ! input parameters
    character(len=*),   parameter :: input_filename = '../../Misc/data/AAAX_16x16/input.txt'
    integer(INT32),     parameter :: SENSOR_DIM1    = 16
    integer(INT32),     parameter :: SENSOR_DIM2    = 16
    
    !character(len=*),   parameter :: input_filename = '../../Misc/data/ABBCBBA_20x20/input.txt'
    !integer(INT32),     parameter :: SENSOR_DIM1    = 20
    !integer(INT32),     parameter :: SENSOR_DIM2    = 20
    
    integer(INT32),     parameter :: N_SENSORS      = SENSOR_DIM1 * SENSOR_DIM2
    
    !========================================================================
    ! runtime parameters
    type param_t
        integer(INT32) :: n_columns
        integer(INT32) :: n_bits_cell
        integer(INT32) :: n_cells_pc
        integer(INT32) :: n_cells
    end type param_t

    integer(INT32),     parameter :: N_COLUMNS      = 500
    integer(INT32),     parameter :: N_BITS_CELL    = 4
    integer(INT8),      parameter :: N_CELLS_PC     = SHIFTL(1, N_BITS_CELL)  ! Number of cells per column
    integer(INT32),     parameter :: N_CELLS        = N_CELLS_PC * N_COLUMNS
    
    ! If true, then during inhibition phase the winning columns are selected as 
    ! the most active columns from the region as a whole. Otherwise, the winning 
    ! columns are selected with respect to their local neighborhoods.
    logical(1),         parameter :: GLOBAL_INHIBITION          = .TRUE.
    
    !========================================================================
    logical(1),         parameter :: DO_POOLING                 = .TRUE.
    
    
    ! the number of new random synapses added to a segment in method replace_synapses
    integer(INT32),     parameter :: NEW_NUMBER_SYNAPSES        = 15

    ! TODO rename to DD_PERMANENCE_INIT
    ! Initial permanence of a new synapse
    real(REAL32),       parameter :: DD_PERMANENCE_INIT         = 0.2
    real(REAL32),       parameter :: SP_PD_PERMANENCE_INIT         = 0.2

    ! TODO rename to DD_PERMANENCE_THRESHOLD
    ! if the permanence value for a synapse is greater than this value, it is said to be connected
    real(REAL32),       parameter :: DD_PERMANENCE_THRESHOLD    = 0.2
    real(REAL32),       parameter :: SP_PD_PERMANENCE_THRESHOLD    = 0.2

    ! Amount by which permanence of synapses are incremented during learning.
    real(REAL32),       parameter :: DD_PERMANENCE_INC       = 0.01
    real(REAL32),       parameter :: PD_PERMANENCE_INC       = 0.01
        
    ! Amount by which permanences of synapses are decremented during learning.
    real(REAL32),       parameter :: DD_PERMANENCE_DEC       = 0.01
    real(REAL32),       parameter :: PD_PERMANENCE_DEC       = 0.01

    ! Misc permanences updates 
    real(REAL32),       parameter :: PD_PERMANENCE_INC_WEAK      = 0.001   ! increment if a column is not often used
    real(REAL32),       parameter :: TP_DD_PREDICTED_SEGMENT_DEC = 0.002   ! decrement if a segment is active incorrectly
    
    ! If the number of active connected synapses on a segment is at least 
    ! this threshold, the segment is said to be active
    integer(INT32),     parameter :: DD_ACTIVATION_THRESHOLD    = 8!number of active synapses needed for TP
    
    ! If the number of synapses active on a segment is at least this 
    ! threshold, it is selected as the best matching cell in a bursting column.
    integer(INT32),     parameter :: MIN_DD_ACTIVATION_THRESHOLD    = 4

    ! number of columns that are active in SP
    integer(INT32),     parameter :: INHIBITION_TOP                = NINT(N_COLUMNS * 0.05)

    ! stimulus_Threshold This is a number specifying the minimum
    ! number of synapses that must be active in order for a column to
    ! turn ON. The purpose of this is to prevent noisy input from
    ! activating columns.
    integer(INT32),     parameter :: STIMULUS_THRESHOLD = 0

    !========================================================================

    integer(INT32),     parameter :: N_DD_SYNAPSES          = 100    ! number of potential synapses per distal dendrite
    integer(INT32),     parameter :: N_PD_SYNAPSES          = 100    ! number of potential synapses per proximal dendrite
    integer(INT32),     parameter :: N_DD_SEGMENTS_MAX      = 2000   ! max number of distal dendrite segments available in a column

    integer(INT32),     parameter :: N_SEGMENT_UPDATES_MAX      = 100
    integer(INT32),     parameter :: SEG_UPDATE_VALID_DURATION  = 5
    type seg_update_t
        integer(INT32)  :: segment              = -1
        integer(INT32)  :: create_date          = -1
        integer(INT8)   :: cell                 = -1
        integer(INT8)   :: n_new_synapes        = 0
        logical(1)      :: sequence_segment     = .FALSE.
        logical(1)      :: dd_synapse_active       (N_DD_SYNAPSES) = .FALSE.
        real(REAL32)    :: dd_synapse_permanence   (N_DD_SYNAPSES) = 0
        integer(INT32)  :: dd_synapse_origin       (N_DD_SYNAPSES) = -1
    end type seg_update_t

    type column_t
        logical(1) :: matching_cells            (N_CELLS_PC)
        logical(1) :: prev_matching_cells       (N_CELLS_PC)
    
        logical(1) :: predicted_inactive_cells  (N_CELLS_PC) ! only needed to prevent allocation in TP
        
        logical(1) :: active_segments           (N_DD_SEGMENTS_MAX)
        logical(1) :: prev_active_segments      (N_DD_SEGMENTS_MAX)
        
        logical(1) :: matching_segments         (N_DD_SEGMENTS_MAX)
        logical(1) :: prev_matching_segments    (N_DD_SEGMENTS_MAX)
        
        logical(1) :: learning_segments         (N_DD_SEGMENTS_MAX) ! only needed to prevent allocation in TP
        logical(1) :: sequence_segments         (N_DD_SEGMENTS_MAX) ! only needed tp:learn_phase_1 
        
        real(REAL32)    :: boost_factor         = 1.0
        
        real(REAL32)    :: pd_synapse_permanence (N_PD_SYNAPSES) = 0
        integer(INT32)  :: pd_synapse_origin     (N_PD_SYNAPSES) = -1
        
        integer(INT32)  :: dd_segment_count        = 0 ! number of segments in use
        integer(INT8)   :: dd_segment_destination  (N_DD_SEGMENTS_MAX) = -1 ! cell to which this segment belongs to
        
        integer(INT32)  :: dd_synapse_count        (N_DD_SEGMENTS_MAX) = 0  ! number of synapses in use
        integer(INT32)  :: dd_synapse_connected_count (N_DD_SEGMENTS_MAX) = 0
        real(REAL32)    :: dd_synapse_permanence   (N_DD_SYNAPSES, N_DD_SEGMENTS_MAX) = 0
        integer(INT32)  :: dd_synapse_origin       (N_DD_SYNAPSES, N_DD_SEGMENTS_MAX) = -1

        type(seg_update_t) seg_updates(N_SEGMENT_UPDATES_MAX)
        
        !TODO rename following 3 to segments instead of synapse
        integer(INT32)  :: dd_synapse_total_activation(N_DD_SEGMENTS_MAX)      = 0
        integer(INT32)  :: dd_synapse_positive_activations(N_DD_SEGMENTS_MAX)  = 0
        integer(INT32)  :: dd_synapse_last_active_iteration(N_DD_SEGMENTS_MAX) = -1
    end type column_t
    type(column_t), allocatable :: columns(:)
    
    ! trace variables
    logical(1), parameter      :: TRACE_ON      = .TRUE.
    integer(INT32), parameter  :: TRACE_COLUMNS = 338
    integer(INT32), parameter  :: TRACE_CELL    = 1

    contains
!========================================================================
    subroutine init_param(n_columns, n_bits_cell, param)
        integer(INT32), intent(in)    :: n_columns, n_bits_cell
        type(param_t),  intent(inout) :: param
    
        param % n_columns   = n_columns
        param % n_bits_cell = n_bits_cell
        param % n_cells_pc  = SHIFTL(1, n_bits_cell)
        param % n_cells     = n_columns * param % n_cells_pc
        
    end subroutine init_param
!========================================================================
end module htm_v3_constants
