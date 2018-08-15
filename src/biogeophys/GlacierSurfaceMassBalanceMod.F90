module GlacierSurfaceMassBalanceMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Computes fluxes that are specific to glaciers
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use clm_varcon     , only : spval, secspday, h2osno_max
  use clm_varpar     , only : nlevgrnd
  use clm_time_manager, only : get_step_size
  use landunit_varcon, only : istice_mec
  use ColumnType     , only : col                
  use LandunitType   , only : lun
  use glc2lndMod     , only : glc2lnd_type
  use WaterstateType , only : waterstate_type
  use WaterfluxType  , only : waterflux_type
  use abortutils     , only : endrun

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  type, public :: glacier_smb_type
     private

     ! ------------------------------------------------------------------------
     ! Public data
     ! ------------------------------------------------------------------------

     real(r8), pointer, public :: qflx_glcice_col(:) ! col net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC; only valid inside the do_smb_c filter
     real(r8), pointer, public :: qflx_glcice_dyn_water_flux_col(:) ! col water flux needed for balance check due to glc_dyn_runoff_routing (mm H2O/s) (positive means addition of water to the system); valid for all columns

     ! ------------------------------------------------------------------------
     ! Private data
     ! ------------------------------------------------------------------------

     real(r8), pointer :: qflx_glcice_frz_col (:)   ! col ice growth (positive definite) (mm H2O/s); only valid inside the do_smb_c filter
     real(r8), pointer :: qflx_glcice_melt_col(:)   ! col ice melt (positive definite) (mm H2O/s); only valid inside the do_smb_c filter


   contains

     ! ------------------------------------------------------------------------
     ! Public routines
     ! ------------------------------------------------------------------------

     procedure, public :: Init                 ! Allocate memory, initialize variables and setup for history and restart
     procedure, public :: Clean                ! Deallocate memory

     ! The science routines need to be separated into a few pieces so they can be
     ! sequenced properly based on what variables they depend on, and what they affect
     procedure, public :: HandleIceMelt        ! compute ice melt in glacier columns, and convert liquid back to ice
     procedure, public :: ComputeSurfaceMassBalance ! compute fluxes other than ice melt
     procedure, public :: AdjustRunoffTerms    ! adjust liquid and ice runoff fluxes due to glacier fluxes

     ! ------------------------------------------------------------------------
     ! Public routines for unit testing
     ! ------------------------------------------------------------------------
     procedure, public :: GetFreeze    ! Return qflx_glcice_frz
     procedure, public :: GetMelt      ! Return qflx_glcice_melt_col
     procedure, public :: CheckNML     ! Check namelist settings
     ! ------------------------------------------------------------------------
     ! Private routines
     ! ------------------------------------------------------------------------

     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, private :: InitReadNML

  end type glacier_smb_type

  !----------------------------------------------------------
  ! glacier_mec control variables: default values (may be overwritten by namelist)
  !----------------------------------------------------------

  ! number of days before one considers the perennially snow-covered point 'land ice'
  ! ONLY PUBLIC FOR UNIT_TESTING!
  integer , public :: glc_snow_persistence_max_days = 7300
  real(r8), public :: glc_snow_min_swe
  logical , public :: glc_smb_include_snowpack      = .true.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Infrastructure routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)
    class(glacier_smb_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*),  intent(in) :: NLFilename
    !-----------------------------------------------------------------------

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)
    glc_snow_min_swe = h2osno_max
    if ( len_trim(NLFilename) > 0 )then
       call this%InitReadNML( NLFilename )
    end if
    call this%CheckNML( )
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    class(glacier_smb_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds

    integer :: begc, endc
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    allocate(this%qflx_glcice_col               (begc:endc)) ; this%qflx_glcice_col                (:) = nan
    allocate(this%qflx_glcice_dyn_water_flux_col(begc:endc)) ; this%qflx_glcice_dyn_water_flux_col (:) = nan
    allocate(this%qflx_glcice_frz_col           (begc:endc)) ; this%qflx_glcice_frz_col            (:) = nan
    allocate(this%qflx_glcice_melt_col          (begc:endc)) ; this%qflx_glcice_melt_col           (:) = nan
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(glacier_smb_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc

    this%qflx_glcice_col(begc:endc) = spval
    call hist_addfld1d (fname='QICE',  units='mm/s',  &
         avgflag='A', long_name='ice growth/melt', &
         ptr_col=this%qflx_glcice_col, l2g_scale_type='ice')

    this%qflx_glcice_frz_col(begc:endc) = spval
    call hist_addfld1d (fname='QICE_FRZ',  units='mm/s',  &
         avgflag='A', long_name='ice growth', &
         ptr_col=this%qflx_glcice_frz_col, l2g_scale_type='ice')

    this%qflx_glcice_melt_col(begc:endc) = spval
    call hist_addfld1d (fname='QICE_MELT',  units='mm/s',  &
         avgflag='A', long_name='ice melt', &
         ptr_col=this%qflx_glcice_melt_col, l2g_scale_type='ice')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    class(glacier_smb_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds

    integer :: c
    !-----------------------------------------------------------------------

    ! Initialize qflx_glcice_dyn_water_flux_col to 0 for all columns because we want this
    ! flux to remain 0 for columns where is is never set, including non-glacier columns.
    !
    ! Other fluxes intentionally remain unset (spval) outside the do_smb filter, so that
    ! they are flagged as missing value outside that filter.
    do c = bounds%begc, bounds%endc
       this%qflx_glcice_dyn_water_flux_col(c) = 0._r8
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine InitReadNML(this, NLFilename)
    ! Read the glacier surface mass balance namelist
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    implicit none
    class(glacier_smb_type), intent(inout) :: this
    character(len=*),  intent(in) :: NLFilename
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'Glacier_SurfaceMassBalance::InitReadNML'
    character(len=*), parameter :: nmlname = 'glacier_smb'

    ! Glacier_mec info
    namelist /glacier_smb/ &
         glc_snow_persistence_max_days, glc_snow_min_swe, glc_smb_include_snowpack

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=glacier_smb, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (glc_snow_persistence_max_days, mpicom)
    call shr_mpi_bcast (glc_snow_min_swe,              mpicom)
    call shr_mpi_bcast (glc_smb_include_snowpack,      mpicom)
  end subroutine InitReadNML

  !-----------------------------------------------------------------------
  subroutine CheckNML(this)
    ! Check the values from the Namelist
    ! !USES:
    use spmdMod        , only : masterproc, mpicom
    use clm_varctl     , only : iulog
    implicit none
    class(glacier_smb_type), intent(inout) :: this

    if ( h2osno_max <= 0.0_r8 )then
       call endrun(msg="ERROR h2osno_max MUST be set before glacier_smb_type Init")
       return
    end if
    if ( glc_snow_persistence_max_days < 0 )then
       call endrun(msg="ERROR glc_snow_persistence_max_days is negative and can not be")
    end if
    if ( glc_snow_min_swe <= 0.0_r8 )then
       call endrun(msg="ERROR glc_snow_min_swe is negative or zero and can not be")
    end if
    if ( glc_snow_min_swe > h2osno_max )then
       call endrun(msg="ERROR glc_snow_min_swe is greater than h2osno_max and can not be")
    end if

    if (masterproc) then
       write(iulog,*) '   glc snow persistence max days      = ', glc_snow_persistence_max_days
       write(iulog,*) '   glc snow min snow water equivalent = ', glc_snow_min_swe
       write(iulog,*) '   glc surf. mass bal. includes snow  = ', glc_smb_include_snowpack
    end if
    !-----------------------------------------------------------------------

  end subroutine CheckNML


  !-----------------------------------------------------------------------
  subroutine Clean(this)
    class(glacier_smb_type), intent(inout) :: this

    !-----------------------------------------------------------------------

    deallocate(this%qflx_glcice_col               )
    deallocate(this%qflx_glcice_dyn_water_flux_col)
    deallocate(this%qflx_glcice_frz_col )
    deallocate(this%qflx_glcice_melt_col)
  end subroutine Clean

  !-----------------------------------------------------------------------
  subroutine GetFreeze( this, bounds, Freeze )
    class(glacier_smb_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), allocatable, intent(out) :: Freeze(:)
    !-----------------------------------------------------------------------
    allocate(Freeze(bounds%begc:bounds%endc))
    Freeze(:) = this%qflx_glcice_frz_col(:)

  end subroutine GetFreeze

  !-----------------------------------------------------------------------
  subroutine GetMelt( this, bounds, Melt )
    class(glacier_smb_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), allocatable, intent(out) :: Melt(:)
    !-----------------------------------------------------------------------
    allocate(Melt(bounds%begc:bounds%endc))
    Melt(:) = this%qflx_glcice_melt_col(:)

  end subroutine GetMelt

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine HandleIceMelt(this, bounds, num_do_smb_c, filter_do_smb_c, &
       waterstate_inst)
    !
    ! !DESCRIPTION:
    ! Compute ice melt in glacier columns, and convert liquid back to ice
    !
    ! Ideally this should be called immediately after ice is melted, so that liquid is
    ! converted back to ice as soon as possible.
    !
    ! NOTE(wjs, 2016-06-29) Currently this is separated from the main ComputeSurfaceMassBalance
    ! routine so that it can be called from the same place in the driver loop where it was
    ! done before the introduction of GlacierSurfaceMassBalanceMod. This was needed to maintain
    ! identical answers, due to the adjustment of h2osoi_ice and h2osoi_liq in this
    ! routine. In principle we should be able to do these adjustments of h2osoi_ice and
    ! h2osoi_liq later in the driver loop: this would just mean that some intervening
    ! science code is operating on the temporarily-thawed state, before the water runs off
    ! and is replaced by ice from below. The main reason to make this change would be to
    ! simplify the driver logic, consolidating calls to this module. On the other hand,
    ! having a period when there is liquid water at the top of the glacier column could
    ! defeat some of the purpose of converting it immediately back to ice (i.e., so that
    ! the surface fluxes are always generated based on an ice-covered surface) - so it
    ! may be best to keep this separate.
    !
    ! !ARGUMENTS:
    class(glacier_smb_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_do_smb_c        ! number of column points in filter_do_smb_c
    integer, intent(in) :: filter_do_smb_c(:)  ! column filter for points where SMB is calculated
    type(waterstate_type), intent(inout) :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: j
    integer :: fc, c, l
    real(r8) :: dtime ! land model time step (sec)

    character(len=*), parameter :: subname = 'HandleIceMelt'
    !-----------------------------------------------------------------------

    associate( &
         qflx_glcice_melt => this%qflx_glcice_melt_col , & ! Output: [real(r8) (:)   ] ice melt (positive definite) (mm H2O/s)
         h2osoi_liq       => waterstate_inst%h2osoi_liq_col      , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
         h2osoi_ice       => waterstate_inst%h2osoi_ice_col        & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
         )

    dtime = get_step_size()

    do fc = 1, num_do_smb_c
       c = filter_do_smb_c(fc)
       qflx_glcice_melt(c) = 0._r8
    end do

    ! Note that, because the following code only operates over the do_smb filter, that
    ! means that the conversion of water back to ice only happens for glacier columns
    ! where we're computing SMB.

    do j = 1, nlevgrnd
       do fc = 1, num_do_smb_c
          c = filter_do_smb_c(fc)
          l = col%landunit(c)

          if (lun%itype(l) == istice_mec) then
             if (h2osoi_liq(c,j) > 0._r8) then   ! ice layer with meltwater
                qflx_glcice_melt(c) = qflx_glcice_melt(c) + h2osoi_liq(c,j)/dtime
                
                ! convert layer back to pure ice by "borrowing" ice from below the column
                h2osoi_ice(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j)
                h2osoi_liq(c,j) = 0._r8
             end if  ! liquid water is present
          end if  ! istice_mec
       end do
    end do

    end associate

  end subroutine HandleIceMelt

  !-----------------------------------------------------------------------
  subroutine ComputeSurfaceMassBalance(this, bounds, num_allc, filter_allc, &
       num_do_smb_c, filter_do_smb_c, glc2lnd_inst, waterstate_inst, waterflux_inst)
    !
    ! !DESCRIPTION:
    ! Compute glacier fluxes other than ice melt.
    !
    ! This sets the public fields qflx_glcice_col and qflx_glcice_dyn_water_flux_col to
    ! their final values.
    !
    ! Should be called after HandleIceMelt, and after waterflux_inst%qflx_snwcp_ice_col is
    ! computed
    !
    ! !ARGUMENTS:
    class(glacier_smb_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_allc        ! number of column points in filter_allc
    integer, intent(in) :: filter_allc(:)  ! column filter for all points
    integer, intent(in) :: num_do_smb_c        ! number of column points in filter_do_smb_c
    integer, intent(in) :: filter_do_smb_c(:)  ! column filter for points where SMB is calculated
    type(glc2lnd_type), intent(in) :: glc2lnd_inst
    type(waterstate_type), intent(in) :: waterstate_inst
    type(waterflux_type), intent(in) :: waterflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, l, g
    real(r8) :: dtime ! land model time step (sec)

    character(len=*), parameter :: subname = 'ComputeSurfaceMassBalance'
    !-----------------------------------------------------------------------

    associate( &
         qflx_glcice                => this%qflx_glcice_col                    , & ! Output: [real(r8) (:)] net flux of new glacial ice (growth - melt) (mm H2O/s)
         qflx_glcice_frz            => this%qflx_glcice_frz_col                , & ! Output: [real(r8) (:)] ice growth (positive definite) (mm H2O/s)
         qflx_glcice_dyn_water_flux => this%qflx_glcice_dyn_water_flux_col     , & ! Output: [real(r8) (:)] water flux needed for balance check due to glc_dyn_runoff_routing (mm H2O/s) (positive means addition of water to the system)
         qflx_glcice_melt           => this%qflx_glcice_melt_col               , & ! Input: [real(r8) (:)] ice melt (positive definite) (mm H2O/s)
         glc_dyn_runoff_routing     => glc2lnd_inst%glc_dyn_runoff_routing_grc , & ! Input:  [real(r8) (:)]  whether we're doing runoff routing appropriate for having a dynamic icesheet
         snow_persistence           => waterstate_inst%snow_persistence_col    , & ! Input: [real(r8) (:)]  counter for length of time snow-covered
         h2osno                     => waterstate_inst%h2osno_col              , & ! Input: [real(r8) (:)] col snow water (mm H2O)
         h2osno_old                 => waterstate_inst%h2osno_old_col          , & ! Input: [real(r8) (:)] col snow mass for previous time step (kg/m2)
         qflx_snwcp_ice             => waterflux_inst%qflx_snwcp_ice_col         & ! Input: [real(r8) (:)]  excess solid h2o due to snow capping (outgoing) (mm H2O /s) [+]
         )
    dtime = get_step_size()

    ! NOTE(wjs, 2016-06-29) The following initialization is done in case the columns
    ! included / excluded in the do_smb_c filter can change mid-run (besides just being
    ! active vs. inactive): If an active column was inside this filter in the previous
    ! timestep, but is no longer inside this filter in this timestep, we want this flux to
    ! be 0 (rather than remaining at its previous value). (Currently, the set of active
    ! columns included in the do_smb filter cannot change mid-run, but the logic is
    ! complex enough that I don't want to assume that that will always remain true.) This
    ! initialization also handles the case where glc_dyn_runoff_routing may change
    ! mid-run, so that a point previously inside that mask no longer is.
    do fc = 1, num_allc
       c = filter_allc(fc)
       qflx_glcice_dyn_water_flux(c) = 0._r8
    end do


    ! Calculate positive surface mass balance to ice sheets, both from already-glaciated
    ! landunits and from non-glaciated landunits (glacial inception)
    do fc = 1, num_do_smb_c
       c = filter_do_smb_c(fc)
       l = col%landunit(c)
       g = col%gridcell(c)
       ! In the following, we convert glc_snow_persistence_max_days to r8 to avoid overflow
       if ( (snow_persistence(c) >= (real(glc_snow_persistence_max_days, r8) * secspday)) &
            .or. lun%itype(l) == istice_mec) then
          qflx_glcice_frz(c) = qflx_snwcp_ice(c)
          ! If SMB for glaciers include snowpack...
          if ( glc_smb_include_snowpack )then
             qflx_glcice_frz(c) = (h2osno(c) - h2osno_old(c)) / dtime + qflx_glcice_frz(c)
          end if
       else
          qflx_glcice_frz(c) = 0._r8
       end if

       qflx_glcice(c) = qflx_glcice_frz(c) - qflx_glcice_melt(c)

       ! For glc_dyn_runoff_routing > 0::
       ! (1) All or part of the excess snow (from snow capping) has been incorporated in
       !     qflx_glcice_frz.  This flux must be included here to complete the water
       !     balance, because it is a sink of water as far as CLM is concerned (this water
       !     will now be owned by CISM).
       ! (2) Meltwater from ice (qflx_glcice_melt) is allowed to run off and is included
       !     in qflx_qrgwl, but the water content of the ice column has not changed
       !     because an equivalent ice mass has been "borrowed" from the base of the
       !     column. So this borrowing is a source of water as far as CLM is concerned.
       !   
       ! For glc_dyn_runoff_routing = 0: Point (2) is the same as for the
       ! glc_dyn_runoff_routing > 0 case: there is a source of water equal to
       ! qflx_glcice_melt. However, in this case, the sink of water is also equal to
       ! qflx_glcice_melt: We have simply replaced some liquid water with an equal amount
       ! of solid ice. Another way to think about this is:
       ! (1) qflx_ice_runoff_snwcp is reduced by an amount equal to qflx_glcice_melt (done
       !     elsewhere in this module). The amount of snow removal is therefore given by
       !     (qflx_ice_runoff_snwcp + qflx_glcice_melt), meaning that there is an
       !     additional sink of water equal to qflx_glcice_melt.
       ! (2) Meltwater from ice (qflx_glcice_melt) is allowed to run off and is included
       !     in qflx_qrgwl, but the water content of the ice column has not changed
       !     because an equivalent ice mass has been "borrowed" from the base of the
       !     column. So this borrowing is a source of water as far as CLM is concerned.
       ! These two corrections cancel out, so nothing is done here.
       qflx_glcice_dyn_water_flux(c) = glc_dyn_runoff_routing(g) * (qflx_glcice_melt(c) - qflx_glcice_frz(c))
    end do

    end associate

  end subroutine ComputeSurfaceMassBalance

  !-----------------------------------------------------------------------
  subroutine AdjustRunoffTerms(this, bounds, num_do_smb_c, filter_do_smb_c, &
       glc2lnd_inst, qflx_qrgwl, qflx_ice_runoff_snwcp)
    !
    ! !DESCRIPTION:
    ! Adjust liquid and ice runoff fluxes due to glacier fluxes
    !
    ! Should be called after ComputeSurfaceMassBalance, and after qflx_qrgwl and
    ! qflx_ice_runoff_snwcp have been given their initial values
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(glacier_smb_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_do_smb_c        ! number of column points in filter_do_smb_c
    integer, intent(in) :: filter_do_smb_c(:)  ! column filter for points where SMB is calculated
    type(glc2lnd_type), intent(in) :: glc2lnd_inst
    real(r8), intent(inout) :: qflx_qrgwl( bounds%begc: ) ! col qflx_surf at glaciers, wetlands, lakes
    real(r8), intent(inout) :: qflx_ice_runoff_snwcp( bounds%begc: ) ! col solid runoff from snow capping (mm H2O /s)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, g

    character(len=*), parameter :: subname = 'AdjustRunoffTerms'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(qflx_qrgwl) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_ice_runoff_snwcp) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         qflx_glcice_frz        => this%qflx_glcice_frz_col                   , & ! Input: [real(r8) (:)] ice growth (positive definite) (mm H2O/s)
         qflx_glcice_melt       => this%qflx_glcice_melt_col                  , & ! Input: [real(r8) (:)] ice melt (positive definite) (mm H2O/s)
         glc_dyn_runoff_routing =>    glc2lnd_inst%glc_dyn_runoff_routing_grc   & ! Input: [real(r8) (:)]  gridcell fraction coupled to dynamic ice sheet
         )

    ! Note that, because the following code only operates over the do_smb filter, that
    ! means that the adjustments here are only applied for glacier columns where we're
    ! computing SMB. This is consistent with the use of the do_smb filter in
    ! HandleIceMelt.

    do fc = 1, num_do_smb_c
       c = filter_do_smb_c(fc)
       g = col%gridcell(c)

       ! Melt is only generated for glacier columns. But it doesn't hurt to do this for
       ! all columns in the do_smb filter: this melt term will be 0 for other columns.
       ! Note: Ice melt is added to the runoff whether or not the column is coupled
       ! to a dynamic glacier model.

       qflx_qrgwl(c) = qflx_qrgwl(c) + qflx_glcice_melt(c)

       ! For the part of the column that is coupled to a dynamic glacier model,
       ! the glacier model handles the fate of capped snow, so we do not want it sent to runoff.
       qflx_ice_runoff_snwcp(c) = qflx_ice_runoff_snwcp(c) - glc_dyn_runoff_routing(g)*qflx_glcice_frz(c)

       ! In places where we are not coupled to a dynamic glacier model, CLM sends all of
       ! the snow capping to the ocean as an ice runoff term. (This is essentially a crude
       ! parameterization of calving, assuming steady state - i.e., all ice gain is
       ! balanced by an ice loss.) But each unit of melt that happens is an indication
       ! that 1 unit of the ice shouldn't have made it to the ocean - but instead melted
       ! before it got there. So we need to correct for that by removing 1 unit of ice
       ! runoff for each unit of melt. Note that, for a given point in space & time, this
       ! can result in negative ice runoff. However, we expect that, in a temporally and
       ! spatially-integrated sense (if we're near equilibrium), this will just serve to
       ! decrease the too-high positive ice runoff.
       !
       ! Another way to think about this is: ice melt removes mass; the snow capping flux
       ! also removes mass. If both the accumulation and melt remove mass, there is a
       ! double-counting. So we need to correct that by: for each unit of ice melt
       ! (resulting in 1 unit of liquid runoff), remove 1 unit of ice runoff. (This is not
       ! an issue for parts of the column coupled to the dynamic glacier model, because
       ! then the snow capping mass is retained in the LND-GLC coupled system.)
       !
       ! The alternative of simply not adding ice melt to the runoff stream where
       ! glc_dyn_runoff_routing = 0 conserves mass, but fails to conserve energy, for a
       ! similar reason: Ice melt in CLM removes energy; also, the ocean's melting of the
       ! snow capping flux removes energy. If both the accumulation and melting remove
       ! energy, there is a double-counting.
       !
       ! Yet another way to think about this is: When ice melted, we let the liquid run
       ! off, and replaced it with new ice from below. But that new ice needed to come
       ! from somewhere to keep the system in water balance. We "request" the new ice from
       ! the ocean by generating a negataive ice runoff equivalent to the amount we have
       ! melted (i.e., equivalent to the amount of new ice that we created from below).
       
       ! As above: Melt is only generated for glacier columns. But it doesn't hurt to do
       ! this for all columns in the do_smb filter: this melt term will be 0 for other
       ! columns.
       
       qflx_ice_runoff_snwcp(c) = qflx_ice_runoff_snwcp(c) - (1.0_r8 - glc_dyn_runoff_routing(g)) * qflx_glcice_melt(c)

       ! Recall that glc_dyn_runoff_routing = min(lfrac, Sg_icemask_coupled_fluxes_l) / lfrac.
       !
       ! Consider a cell with lfrac = 0.8 and Sg_icemask_coupled_fluxes_l = 0.4. (For
       ! instance, the cell might have half its land coverage in Greenland and the other
       ! half in Ellemere.)  Given qflx_ice_runoff_snwcp(c) = 1 m/yr, half the flux (0.5
       ! m/yr) will be sent to the runoff model, where it will be multiplied by lfrac to
       ! give a net flux of 0.4 m/yr times the cell area.
       !
       ! The full SMB of 1 m/yr will be sent to the coupler's prep_glc_mod, but it will be
       ! weighted by 0.4 when integrating over the whole ice sheet. So a net flux of 0.4
       ! m/yr will also be applied to the ice sheet model.  The total flux of 0.8 m/yr,
       ! split evenly between runoff and ice sheet, is what we want.

    end do

    end associate

  end subroutine AdjustRunoffTerms

end module GlacierSurfaceMassBalanceMod
