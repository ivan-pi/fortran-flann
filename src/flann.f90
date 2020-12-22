module flann

  use iso_c_binding, only: c_ptr, c_float, c_int
  use flann_c
  use flann_params

  implicit none
  private

  integer, parameter :: wp = c_float

  !
  ! Exposed constant from flann_c
  !
  !
  ! flann_centers_init_t
  public :: FLANN_CENTERS_RANDOM
  public :: FLANN_CENTERS_GONZALES
  public :: FLANN_CENTERS_KMEANSPP
  public :: FLANN_CENTERS_GROUPWISE

  !
  ! flann_log_level_t
  public :: FLANN_LOG_NONE
  public :: FLANN_LOG_FATAL
  public :: FLANN_LOG_ERROR
  public :: FLANN_LOG_WARN
  public :: FLANN_LOG_INFO
  public :: FLANN_LOG_DEBUG

  !
  ! flann_distance_t
  public :: FLANN_DIST_EUCLIDEAN       
  public :: FLANN_DIST_L2              
  public :: FLANN_DIST_MANHATTAN       
  public :: FLANN_DIST_L1              
  public :: FLANN_DIST_MINKOWSKI       
  public :: FLANN_DIST_MAX             
  public :: FLANN_DIST_HIST_INTERSECT  
  public :: FLANN_DIST_HELLINGER       
  public :: FLANN_DIST_CHI_SQUARE      
  public :: FLANN_DIST_KULLBACK_LEIBLER
  public :: FLANN_DIST_HAMMING         
  public :: FLANN_DIST_HAMMING_LUT     
  public :: FLANN_DIST_HAMMING_POPCNT  
  public :: FLANN_DIST_L2_SIMPLE       

  !
  ! Exposed entities from `flann_params` module.
  !
  public :: flann_log_verbosity
    !! Set the log verbosity level.
  public :: search_params
  public :: LinearIndexParams
  public :: KDTreeIndexParams
  public :: KMeansIndexParams
  public :: CompositeIndexParams
  public :: KDTreeSingleIndexParams
  public :: HierarchicalClusteringIndexParams
  public :: LshIndexParams
  public :: AutotunedIndexParams

  !
  ! Public entities from this module
  !
  public :: flann_index
    !! The FLANN nearest neighbor type.
  public :: flann_hierarchical_clustering
    !! Clusters the given points by constructing a hierarchical k-means tree
    !! and choosing a cut in the tree that minimized the clusters' variance.


  !> The FLANN nearest neighbor index type.
  !
  ! This derived type is used to abstract different types of
  ! nearest neighbor search indexes.
  type flann_index
    type(c_ptr) :: index_id = c_null_ptr
    type(FLANNParameters) :: params
    real(wp), pointer :: points(:,:) => null()
    integer(c_int) :: rows, cols
  contains
    private
    procedure :: build_index_v1
    procedure :: build_index_v2
    generic, public :: build_index => build_index_v1, build_index_v2
      !! Build the nearest neighbor index. There are two versions of
      !! this method, one that uses the points provided as argument and
      !! one that uses the the points provided when the object was
      !! constructed.
    procedure, public :: knn_search => index_knn_search
      !! Performs a K-nearest neighbor search for a set of query points.
    procedure, public :: radius_search => index_radius_search
      !! Performs a radius nearest neighbor search for a set of query points.
    procedure, public :: save => index_save
      !! Saves the index to a file.
    procedure, public :: veclen => index_veclen
      !! Returns number of features in this index.
    procedure, public :: size => index_size
      !! Returns the dimensionality of the points in this index.
    procedure, public :: get_type => index_get_type
      !! Returns the index type (kdtree, kmeans, ...).
    procedure, public :: used_memory => index_used_memory
      !! Returns the amount of memory (in bytes) used by the index.
    procedure, public :: get_parameters => index_get_parameters
      !! Returns the Index Parameters.
    final :: index_free
  end type

  interface flann_index
    !! Overload the structure constructor.
    module procedure index_new_int
    module procedure index_new_char
    module procedure index_new_int_points
    module procedure index_new_char_points
  end interface

contains

  function distance_helper(distance) result(d)
    character(len=*), intent(in) :: distance
    integer(c_int) :: d

    character(len=:), allocatable :: upper
    upper = uppercase(trim(distance))

    select case(upper)
    case('EUCLIDEAN','L2')
      d = FLANN_DIST_EUCLIDEAN
    case('MANHATTAN','L1')
      d = FLANN_DIST_MANHATTAN
    case('MINKOWSKI')
      d = FLANN_DIST_MINKOWSKI
    case('MAX_DIST')
      d = FLANN_DIST_MAX
    case('HIK','HIST_INTERSECT')
      d = FLANN_DIST_HIST_INTERSECT
    case('HELLINGER')
      d = FLANN_DIST_HELLINGER
    case('CS','CHI_SQUARE')
      d = FLANN_DIST_CHI_SQUARE
    case('KL','KULLBACK_LEIBLER')
      d = FLANN_DIST_KULLBACK_LEIBLER
    ! case('HAMMING')
    !   d = FLANN_DIST_HAMMING
    ! case('HAMMING_LUT')
    !   d = FLANN_DIST_HAMMING_LUT
    ! case('HAMMING_POPCNT')
    !   d = FLANN_DIST_HAMMING_POPCNT
    case('L2_SIMPLE')
      d = FLANN_DIST_L2_SIMPLE
    case default
      ! Unsupported flann_distance_t (an error will be raised on the C interface level.)
      d = 999 
    end select    

  contains

    function uppercase(str) result(upper_str)
      character(len=*), intent(in) :: str
      character(len=len(str)) :: upper_str
      integer :: i, ic, k
      character(len=26), parameter :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(len=26), parameter :: lower_case = 'abcdefghijklmnopqrstuvwxyz'

      do i = 1, len(str)
        ic = iachar(str(i:i))
        if (ic >= iachar('a') .and. ic <= iachar('z')) then
          ! lower case
          k = index(lower_case,str(i:i))
          upper_str(i:i) = upper_case(k:k)
        else
          ! upper case or non-alphabetic
          upper_str(i:i) = str(i:i)
        end if
      end do
    end function

  end function

  ! function Index_empty_char(params,distance) result(this)
  !   type(FLANNParameters), intent(in) :: params
  !   character(len=*), intent(in) :: distance
  !   type(flann_index) :: this



  ! end function

  function index_new_int(params,distance,order) result(this)
    type(FLANNParameters), intent(in) :: params
    integer(c_int), intent(in) :: distance
    integer(c_int), intent(in), optional :: order
    type(flann_index) :: this
    integer(c_int) :: order_

    this%params = params

    order_ = 2 ! Euclidean distance
    if (present(order)) order_ = order

    call flann_set_distance_type(distance,order_)
  end function

  function index_new_char(params,distance,order) result(this)
    type(FLANNParameters), intent(in) :: params
    character(len=*), intent(in) :: distance
    integer(c_int), intent(in), optional :: order
    type(flann_index) :: this
    this = index_new_int(params,distance_helper(distance),order)
  end function

  function index_new_int_points(points,params,distance,order) result(this)
    real(wp), intent(in), target :: points(:,:)
      !! Array containing the points (features) that should be indexed,
      !! stored in a column-major order (one point/feature per column).
      !! The shape of the matrix is *(dimensionality)* × *(number of features)*.
    type(FLANNParameters), intent(in) :: params
      !! Structure containing the index parameters.
    integer(c_int), intent(in) :: distance
    integer(c_int), intent(in), optional :: order

    type(flann_index) :: this
    integer(c_int) :: order_
    intrinsic :: size

    this%points => points
    this%rows = size(this%points,1)
    this%cols = size(this%points,2)

    this%params = params

    order_ = 2 ! Euclidean distance
    if (present(order)) order_ = order
    call flann_set_distance_type(distance,order_)

  end function

  function index_new_char_points(points,params,distance,order) result(this)
    real(wp), intent(in), target :: points(:,:)
      !! Array containing the points (features) that should be indexed,
      !! stored in a column-major order (one point/feature per column).
      !! The shape of the matrix is *(dimensionality)* × *(number of features)*.
    type(FLANNParameters), intent(in) :: params
      !! Structure containing the index parameters.
    character(len=*), intent(in) :: distance
    integer(c_int), intent(in), optional :: order

    type(flann_index) :: this

    this = index_new_int_points(points,params,distance_helper(distance),order)
  end function

  subroutine build_index_v1(self,speedup)
    class(flann_index), intent(inout) :: self
    real(c_float), intent(out), optional :: speedup

    real(c_float) :: speedup_

    associate(points => self%points, &
              rows => self%rows, &
              cols => self%cols, &
              params => self%params)
      self%index_id = flann_build_index_float(points,cols,rows,speedup_,params)
    end associate
    if (present(speedup)) speedup = speedup_
  end subroutine

  subroutine build_index_v2(self,points,speedup)
    class(flann_index), intent(inout) :: self
    real(wp), intent(in), target :: points(:,:)
    real(c_float), intent(out), optional :: speedup

    real(c_float) :: speedup_
    intrinsic :: size

    self%points => points
    self%rows = size(self%points,1)
    self%cols = size(self%points,2)    

    associate(points => self%points, &
              rows => self%rows, &
              cols => self%cols, &
              params => self%params)
      self%index_id = flann_build_index_float(points,cols,rows,speedup_,params)
    end associate
    if (present(speedup)) speedup = speedup_
  end subroutine

  subroutine index_knn_search(self,queries,indices,dists,knn,params,stat)
    class(flann_index), intent(inout) :: self
    real(wp), intent(in) :: queries(:,:)
      !! Array containing the query points. 
      !! Size of the array is *(dimensionality)* × *(number of queries)* 
    integer(c_int), intent(inout) :: indices(:,:)
      !! Array that will contain the distances to the K-nearest neighbors found
      !! (size should be at least *(knn)* × *(number of queries)*).
    real(wp), intent(inout) :: dists(:,:)
      !! Array that will contain the distances to the K-nearest neighbors found
      !! (size should be at least *(knn)* × *(number of queries)*).
      !! Note: For Euclidean distances, FLANN computes squared distances, so the 
      !! values returned here are squared distances as well.
    integer(c_int), intent(in) :: knn
      !! Number of nearest neighbors to search for.
    type(search_params), intent(in) :: params
      !! Search parameters. Structure containing parameters used during search.
    integer(c_int), intent(out), optional :: stat

    integer(c_int) :: stat_, qcols

    call copy_search_params(self%params,params)
    
    qcols = size(queries,2)

    if (size(indices,1) /= knn) error stop "Error in knnSearch 1"
    if (size(indices,2) /= qcols) error stop "Error in knnSearch 2"
    if (any(shape(indices) /= shape(dists))) error stop "Error in knnSearch 3"
    ! assert(size(indices,1) == knn)
    ! assert(size(indices,2) == qcols)
    ! assert(all(shape(indices) == shape(dists)))

    print *, c_associated(self%index_id)
    print *, "veclen", self%veclen(), "size ", self%size(), "qdims ", size(queries,1)
    

    print *,associated(self%points), shape(self%points)

    stat_ = flann_find_nearest_neighbors_index_float(self%index_id, &
      queries,qcols,indices,dists,knn,self%params)

    if (present(stat)) then
      stat = stat_
    else
      if (stat_ /= 0) then
        error stop "Error in knnSearch"
      end if
    end if

  end subroutine

  subroutine index_radius_search(self,query,radius,indices,dists,nfound,params,stat)
    class(flann_index), intent(inout) :: self
    real(wp), intent(in) :: query(:)
      !! A single query point. 
    real(wp), intent(in) :: radius
      !! Search radius (squared radius for Euclidean metric).
    integer(c_int), intent(out) :: indices(:)
    real(wp), intent(out) :: dists(size(indices))
    integer(c_int), intent(out) :: nfound
    type(search_params), intent(in) :: params
    integer(c_int), intent(out), optional :: stat

    integer(c_int) :: stat_, max_nn

    ! assert(size(indices) == size(dists))
    ! assert(size(queries) = self%rows)

    max_nn = size(indices)
    call copy_search_params(self%params,params)

    associate(index_id => self%index_id, &
              params => self%params)
      stat_ = flann_radius_search(index_id,query,indices,dists,max_nn,radius,params)
    end associate
    nfound = stat_
    if (present(stat)) then
      stat = stat_
    else
      if (stat < 0) then
        error stop "Error during radiusSearch."
      end if
    end if
  end subroutine

  !> Saves the index to a file.
  !
  subroutine index_save(self,filename,stat)
    use iso_c_binding, only: c_null_char
    class(flann_index), intent(in) :: self
    character(len=*), intent(in) :: filename
      !! The file to save the index to
    integer(c_int), intent(out), optional :: stat
      !! Status flag. Zero on success, negative value on error.
    integer(c_int) :: stat_
    
    ! int flann_save_index(flann_index_t index_id,char* filename);
    stat_ = flann_save_index_float(self%index_id,filename//c_null_char)
    if(present(stat))  then
      stat_ = stat
    else
      if (stat_ /= 0) then
        error stop "Error saving file."
      end if
    end if

  end subroutine

  subroutine index_free(self)
    type(flann_index), intent(inout) :: self

    integer(c_int) :: stat_

    if(associated(self%points)) nullify(self%points)

    if (c_associated(self%index_id)) then
      stat_ = flann_free_index_float(self%index_id, self%params)
      if (stat_ /= 0) then
        write(*,*) "An error occured while freeing the index_id."
        write(*,*) "Stat = ", stat_
      end if
    end if
  end subroutine

  !> Returns number of features (points) in this index.
  !
  integer(c_int) function index_veclen(self) result(veclen)
    class(flann_index), intent(in) :: self
    veclen = self%cols
  end function

  !> Returns the dimensionality of the points in this index.
  !
  integer(c_int) function index_size(self) result(size)
    class(flann_index), intent(in) :: self
    size = self%rows
  end function

  !> Returns the index type (kdtree, kmeans, ...).
  !
  integer(c_int) function index_get_type(self) result(type)
    class(flann_index), intent(in) :: self
    type = self%params%algorithm
  end function

  !> Returns the amount of memory (in bytes) used by the index.
  !
  integer(c_int) function index_used_memory(self) result(used_memory)
    class(flann_index), intent(in) :: self

    ! Retrieve the amount of memory on the C side.
    used_memory = flann_used_memory_float(self%index_id)

    ! Add the size of the meta-data stored on Fortran side
    ! TODO: Verify if storage_size returns a multiple of 8.
    used_memory = used_memory + storage_size(self)/8
  end function

  !> Returns the Index Parameters.
  !
  type(FLANNParameters) function index_get_parameters(self) result(params)
    class(flann_index), intent(in) :: self
    params = self%params
  end function

  !> Clusters the given points by constructing a hierarchical k-means tree and 
  ! and choosing a cut in the tree that minimizes the clusters' variance.
  !
  subroutine flann_hierarchical_clustering(points,clusters,centers,params,distance,order,stat)
    real(wp), intent(in) :: points(:,:)
      !! The points to be clustered.
    integer, intent(inout) :: clusters
      !! On input contains the number of clusters desired. On output contains
      !! the number of clusters computed. If an error occured a value < 0 is returned.
    real(wp), intent(out) :: centers(size(points,1),clusters)
      !! The centers of the clusters obtained. The number of rows
      !! in this array should match the number of clusters desired
      !! (`size(centers,1) == clusters`). However, because of the way
      !! the cut in the hierarchical tree is choosen, the number of clusters
      !! computer will be the highest number of the form `(branching - 1)*k + 1`
      !! that's lower than the number of clusters desired, where `branching` is
      !! is the tree's branching factor.
    type(FLANNParameters), intent(in) :: params
      !! Parameters used in the construction of the hierarchical k-means tree
      !! using KMeansIndexParams().
    integer(c_int), intent(in) :: distance
    integer(c_int), intent(in), optional :: order
    integer(c_int), intent(out), optional :: stat

    integer(c_int) :: rows, cols, order_

    rows = size(points,1)
    cols = size(points,2)

    order_ = 2
    if (present(order)) order_ = order
    call flann_set_distance_type(distance,order_)

    clusters = flann_compute_cluster_centers(points,rows,cols,clusters,centers,params)

    if (present(stat)) then
      if (clusters >= 0) then
        stat = 0
      else
        stat = clusters
      end if
    else
      if (clusters < 0) then
        error stop "Error computing clusters."
      end if
    end if
  end subroutine

end module
