module flann_params

  use iso_c_binding, only: c_int, c_float
  use flann_c
  
  implicit none
  private

  public :: log_verbosity

  public :: search_params, copy_search_params

  public :: LinearIndexParams
  public :: KDTreeIndexParams
  public :: KMeansIndexParams
  public :: CompositeIndexParams
  public :: KDTreeSingleIndexParams
  public :: HierarchicalClusteringIndexParams
  public :: LshIndexParams
  public :: AutotunedIndexParams

  type :: search_params
    integer(c_int) :: checks = 32
      !! Specified the maximum leafs to visit when searching for neighbors. 
      !! A higher value for this parameter would give better search precision, 
      !! but also take more time. For all leafs to be checked use the value 
      !! `FLANN_CHECKS_UNLIMITED`. If automatic configuration was used when
      !! the index was created, the number of checks required to achieve the
      !! specified precision was also computed, to use that value specify
      !! `FLANN_CHECKS_AUTOTUNED`.
    real(c_float) :: eps = 0.0_c_float
      !! Search for eps-approximate neighbors (only used by KDTreeSingleIndex).
    logical :: sorted = .true.
      !! Used only by radius search, specifies if the neighbors
      !! returned should be sorted by distance. 
    integer(c_int) :: max_neighbors = -1
      !! Specifies the maximum number of neighbors radius
      !! search should return (default: -1 = unlimited). Only 
      !! used for radius search.
    integer(c_int) :: cores = 0
      !! How many cores to assign to the search (specify 0 for
      !! automatic core selection).
  end type

  ! interface search_params
  !   module procedure search_params_constructor
  ! end interface

  integer(c_int), protected :: GLOBAL_log_level = FLANN_LOG_NONE

contains

  subroutine log_verbosity(level)
    use flann_c, only: flann_log_verbosity
    integer, intent(in) :: level

    ! Store a copy of the log verbosity in order to
    ! be able to initialize the FLANN Parameters struct
    GLOBAL_log_level = level
    
    ! Set the verbosity on the C side
    call flann_log_verbosity(level)
  end subroutine

  ! function search_params_constructor(checks,eps,sorted) result(this)
  !   integer(c_int), intent(in), optional :: checks
  !   real(c_float), intent(in), optional :: eps
  !   logical, intent(in), optional :: sorted
  !   type(search_params) :: this
  !   if (present(checks)) this%checks = checks
  !   if (present(checks)) this%eps = eps
  !   if (present(sorted)) this%sorted = sorted
  ! end function

  subroutine copy_search_params(params,params_in)
    type(FLANNParameters), intent(inout) :: params
    type(search_params), intent(in) :: params_in

    params%checks = params_in%checks
    params%eps = params_in%eps
    ! if (params_in%sorted) then
    !   params%sorted = 1_c_int
    ! else
    !   params%sorted = 0_c_int
    ! end if
    params%max_neighbors = params_in%max_neighbors
    params%cores = params_in%cores
  end subroutine

  function LinearIndexParams() result(this)
    type(FLANNParameters) :: this
    this%algorithm = FLANN_INDEX_LINEAR
    this%log_level = GLOBAL_log_level
  end function

  function KDTreeIndexParams(trees) result(this)
    integer(c_int), intent(in), optional :: trees
      !! The number of parallel kd-trees to use. Good values are in the 
      !! range from 1 to 16.
    type(FLANNParameters) :: this
    this%algorithm = FLANN_INDEX_KDTREE
    if (present(trees)) this%trees = trees
    this%log_level = GLOBAL_log_level
  end function

  function KMeansIndexParams(branching,iterations,centers_init,cb_index) result(this)
    integer(c_int), intent(in), optional :: branching
      !! The branching factor to use for the hierarchical k-means tree.
    integer(c_int), intent(in), optional :: iterations
      !! The maximum number of iterations to use in the k-means clustering
      !! stage when building the k-means tree. If a value of -1 is used here,
      !! it means that the k-means tree clustering should be iterated until
      !! convergence.
    integer(c_int), intent(in), optional :: centers_init
      !! The algorithm to use when performing a k-means clustering step. The 
      !! possible values are:
      !!
    real(c_float), intent(in), optional :: cb_index
      !! This parameter (cluster boundary index) influences the way
      !! exploration is performed in the hierarchical k-means tree. 
      !! When cb_index is zero the next k-means domain to be explores is 
      !! choosen to be the one with the closest center. A value greater then
      !! zero also takes into account the size of the domain.

    type(FLANNParameters) :: this

    this%algorithm = FLANN_INDEX_KMEANS
    if (present(branching)) this%branching = branching
    if (present(iterations)) this%iterations = iterations
    if (present(centers_init)) this%centers_init = centers_init
    if (present(cb_index)) this%cb_index = cb_index
    this%log_level = GLOBAL_log_level
  end function

  function CompositeIndexParams(trees,branching,iterations,centers_init,cb_index) result(this)
    integer(c_int), intent(in), optional :: trees
    integer(c_int), intent(in), optional :: branching
    integer(c_int), intent(in), optional :: iterations
    integer(c_int), intent(in), optional :: centers_init
    real(c_float), intent(in), optional :: cb_index
    type(FLANNParameters) :: this

    this%algorithm = FLANN_INDEX_COMPOSITE
    if (present(trees)) this%trees = trees
    if (present(branching)) this%branching = branching
    if (present(iterations)) this%iterations = iterations
    if (present(centers_init)) this%centers_init = centers_init
    if (present(cb_index)) this%cb_index = cb_index
    this%log_level = GLOBAL_log_level
  end function

  function KDTreeSingleIndexParams(leaf_max_size) result(this)
    integer(c_int), intent(in), optional :: leaf_max_size
    type(FLANNParameters) :: this
    this%algorithm = FLANN_INDEX_KDTREE_SINGLE
    if (present(leaf_max_size)) this%leaf_max_size = leaf_max_size
    this%log_level = GLOBAL_log_level
  end function

  function HierarchicalClusteringIndexParams(branching,centers_init,trees,leaf_max_size) result(this)
    integer(c_int), intent(in), optional :: branching
    integer(c_int), intent(in), optional :: centers_init
    integer(c_int), intent(in), optional :: trees
    integer(c_int), intent(in), optional :: leaf_max_size
    type(FLANNParameters) :: this
    this%algorithm = FLANN_INDEX_HIERARCHICAL
    if (present(branching)) this%branching = branching
    if (present(centers_init)) this%centers_init = centers_init
    if (present(trees)) this%trees = trees
    if (present(leaf_max_size)) this%leaf_max_size = leaf_max_size
    this%log_level = GLOBAL_log_level
  end function

  function LshIndexParams(table_number,key_size,multi_probe_level) result(this)
    integer(c_int), intent(in), optional :: table_number
    integer(c_int), intent(in), optional :: key_size
    integer(c_int), intent(in), optional :: multi_probe_level
    type(FLANNParameters) :: this
    this%algorithm = FLANN_INDEX_LSH
    if (present(table_number)) this%table_number_ = table_number
    if (present(key_size)) this%key_size_ = key_size
    if (present(multi_probe_level)) this%multi_probe_level_ = multi_probe_level
    this%log_level = GLOBAL_log_level
  end function

  function AutotunedIndexParams(target_precision,build_weight,memory_weight,sample_fraction) result(this)
    real(c_float), intent(in), optional :: target_precision 
    real(c_float), intent(in), optional :: build_weight
    real(c_float), intent(in), optional :: memory_weight
    real(c_float), intent(in), optional :: sample_fraction
    type(FLANNParameters) :: this
    this%algorithm = FLANN_INDEX_AUTOTUNED
    if (present(target_precision)) this%target_precision = target_precision
    if (present(build_weight)) this%build_weight = build_weight
    if (present(memory_weight)) this%memory_weight = memory_weight
    if (present(sample_fraction)) this%sample_fraction = sample_fraction
    this%log_level = GLOBAL_log_level
  end function

end module