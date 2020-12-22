module flann_c

  use iso_c_binding
  implicit none
  
  public

  !
  ! Nearest neighbour index algorithms
  ! flann_algorithm_t
  integer(c_int), parameter :: FLANN_INDEX_LINEAR          = 0
  integer(c_int), parameter :: FLANN_INDEX_KDTREE          = 1
  integer(c_int), parameter :: FLANN_INDEX_KMEANS          = 2
  integer(c_int), parameter :: FLANN_INDEX_COMPOSITE       = 3
  integer(c_int), parameter :: FLANN_INDEX_KDTREE_SINGLE   = 4
  integer(c_int), parameter :: FLANN_INDEX_HIERARCHICAL    = 5
  integer(c_int), parameter :: FLANN_INDEX_LSH             = 6
  integer(c_int), parameter :: FLANN_INDEX_KDTREE_CUDA     = 7
  integer(c_int), parameter :: FLANN_INDEX_SAVED           = 254
  integer(c_int), parameter :: FLANN_INDEX_AUTOTUNED       = 255

  !
  ! flann_centers_init_t
  integer(c_int), parameter :: FLANN_CENTERS_RANDOM = 0
  integer(c_int), parameter :: FLANN_CENTERS_GONZALES = 1
  integer(c_int), parameter :: FLANN_CENTERS_KMEANSPP = 2
  integer(c_int), parameter :: FLANN_CENTERS_GROUPWISE = 3

  !
  ! flann_log_level_t
  integer(c_int), parameter :: FLANN_LOG_NONE = 0
  integer(c_int), parameter :: FLANN_LOG_FATAL = 1
  integer(c_int), parameter :: FLANN_LOG_ERROR = 2
  integer(c_int), parameter :: FLANN_LOG_WARN = 3
  integer(c_int), parameter :: FLANN_LOG_INFO = 4
  integer(c_int), parameter :: FLANN_LOG_DEBUG = 5

  !
  ! flann_distance_t
  integer(c_int), parameter :: FLANN_DIST_EUCLIDEAN            = 1
  integer(c_int), parameter :: FLANN_DIST_L2                   = 1
  integer(c_int), parameter :: FLANN_DIST_MANHATTAN            = 2
  integer(c_int), parameter :: FLANN_DIST_L1                   = 2
  integer(c_int), parameter :: FLANN_DIST_MINKOWSKI            = 3
  integer(c_int), parameter :: FLANN_DIST_MAX                  = 4
  integer(c_int), parameter :: FLANN_DIST_HIST_INTERSECT       = 5
  integer(c_int), parameter :: FLANN_DIST_HELLINGER            = 6
  integer(c_int), parameter :: FLANN_DIST_CHI_SQUARE           = 7
  integer(c_int), parameter :: FLANN_DIST_KULLBACK_LEIBLER     = 8
  integer(c_int), parameter :: FLANN_DIST_HAMMING              = 9
  integer(c_int), parameter :: FLANN_DIST_HAMMING_LUT          = 10
  integer(c_int), parameter :: FLANN_DIST_HAMMING_POPCNT       = 11
  integer(c_int), parameter :: FLANN_DIST_L2_SIMPLE            = 12

  !
  ! flann_datatype_t
  integer(c_int), parameter :: FLANN_NONE      = -1
  integer(c_int), parameter :: FLANN_INT8      = 0
  integer(c_int), parameter :: FLANN_INT16     = 1
  integer(c_int), parameter :: FLANN_INT32     = 2
  integer(c_int), parameter :: FLANN_INT64     = 3
  integer(c_int), parameter :: FLANN_UINT8     = 4
  integer(c_int), parameter :: FLANN_UINT16    = 5
  integer(c_int), parameter :: FLANN_UINT32    = 6
  integer(c_int), parameter :: FLANN_UINT64    = 7
  integer(c_int), parameter :: FLANN_FLOAT32   = 8
  integer(c_int), parameter :: FLANN_FLOAT64   = 9

  !
  ! flann_checks_t
  integer(c_int), parameter :: FLANN_CHECKS_UNLIMITED = -1
  integer(c_int), parameter :: FLANN_CHECKS_AUTOTUNED = -2


  type, bind(c) :: FLANNParameters

    integer(c_int) :: algorithm = FLANN_INDEX_KDTREE ! the algorithm to use

    ! search time parameters
    integer(c_int) :: checks = 32           ! how many leafs (features) to check in one search
    real(c_float) :: eps = 0.0_c_float      ! eps parameter for eps-knn search
    integer(c_int) :: sorted = 0            ! indicates if results returned by radius search should be sorted or not
    integer(c_int) :: max_neighbors = -1    ! limits the maximum number of neighbors should be returned by radius search
    integer(c_int) :: cores = 0             ! number of paralel cores to use for searching

    ! kdtree index parameters
    integer(c_int) :: trees = 4             ! number of randomized trees to use (for kdtree)
    integer(c_int) :: leaf_max_size = 4

    ! kmeans index parameters
    integer(c_int) :: branching = 32        ! branching factor (for kmeans tree)
    integer(c_int) :: iterations = 11       ! max iterations to perform in one kmeans cluetering (kmeans tree)
    integer(c_int) :: centers_init = FLANN_CENTERS_RANDOM   ! algorithm used for picking the initial cluster centers for kmeans tree
    real(c_float) :: cb_index = 0.2_c_float ! cluster boundary index. Used when searching the kmeans tree

    ! autotuned index parameters
    real(c_float) :: target_precision = 0.9_c_float     ! precision desired (used for autotuning, -1 otherwise)
    real(c_float) :: build_weight = 0.01_c_float        ! build tree time weighting factor
    real(c_float) :: memory_weight = 0                  ! index memory weigthing factor
    real(c_float) :: sample_fraction = 0.1_c_float      ! what fraction of the dataset to use for autotuning

    ! LSH parameters
    integer(c_int) :: table_number_ = 0         ! The number of hash tables to use
    integer(c_int) :: key_size_ = 0             ! The length of the key in the hash tables
    integer(c_int) :: multi_probe_level_ = 0    ! Number of levels to use in multi-probe LSH, 0 for standard LSH

    ! other parameters
    integer(c_int) :: log_level = FLANN_LOG_NONE    ! determines the verbosity of each flann function
    integer(c_long) :: random_seed = 0              ! random seed to use
  end type


  interface
    ! 
    ! Sets the log level used for all flann functions (unless
    ! specified in FLANNParameters for each call
    !
    ! Params:
    !  level = verbosity level
    ! 
    subroutine flann_log_verbosity(level) bind(c,name="flann_log_verbosity")
        import c_int
        integer(c_int), value :: level
    end subroutine

    ! 
    ! Sets the distance type to use throughout FLANN.
    ! If distance type specified is MINKOWSKI, the second argument
    ! specifies which order the minkowski distance should have.
    ! 
    subroutine flann_set_distance_type(distance_type,order) bind(c,name="flann_set_distance_type")
        import c_int
        integer(c_int), value :: distance_type
        integer(c_int), value :: order
    end subroutine

    ! 
    ! Gets the distance type in use throughout FLANN.
    ! 
    integer(c_int) function flann_get_distance_type() bind(c,name="flann_get_distance_type")
        import c_int
    end function

    ! 
    ! Gets the distance order in use throughout FLANN (only applicable if minkowski distance
    ! is in use).
    ! 
    integer(c_int) function flann_get_distance_order() bind(c,name="flann_get_distance_order")
        import c_int
    end function

    !
    ! Builds and returns an index. It uses autotuning if the target_precision field of index_params
    ! is between 0 and 1, or the parameters specified if it's -1.
    !
    ! Params:
    !  dataset = pointer to a data set stored in row major order
    !  rows = number of rows (features) in the dataset
    !  cols = number of columns in the dataset (feature dimensionality)
    !  speedup = speedup over linear search, estimated if using autotuning, output parameter
    !  index_params = index related parameters
    !  flann_params = generic flann parameters
    !
    ! Returns: the newly created index or a number <0 for error
    !
    type(c_ptr) function flann_build_index(dataset,rows,cols,speedup,flann_params) bind(c,name="flann_build_index")
        import c_float, c_int, c_ptr, FLANNParameters
        integer(c_int), intent(in), value :: rows
        integer(c_int), intent(in), value :: cols
        real(c_float), intent(in) :: dataset(cols,rows)
        real(c_float), intent(out) :: speedup
        type(FLANNParameters), intent(in) :: flann_params
    end function
    type(c_ptr) function flann_build_index_float(dataset,rows,cols,speedup,flann_params) bind(c,name="flann_build_index_float")
        import c_float, c_int, c_ptr, FLANNParameters
        integer(c_int), intent(in), value :: rows
        integer(c_int), intent(in), value :: cols
        real(c_float), intent(in) :: dataset(cols,rows)
        real(c_float), intent(out) :: speedup
        type(FLANNParameters), intent(in) :: flann_params
    end function
    type(c_ptr) function flann_build_index_double(dataset,rows,cols,speedup,flann_params) bind(c,name="flann_build_index_float")
        import c_float, c_int, c_ptr, c_double, FLANNParameters
        integer(c_int), intent(in), value :: rows
        integer(c_int), intent(in), value :: cols
        real(c_double), intent(in) :: dataset(cols,rows)
        real(c_float), intent(out) :: speedup
        type(FLANNParameters), intent(in) :: flann_params
    end function
    type(c_ptr) function flann_build_index_int(dataset,rows,cols,speedup,flann_params) bind(c,name="flann_build_index_float")
        import c_float, c_int, c_ptr, FLANNParameters
        integer(c_int), intent(in), value :: rows
        integer(c_int), intent(in), value :: cols
        integer(c_int), intent(in) :: dataset(cols,rows)
        real(c_float), intent(out) :: speedup
        type(FLANNParameters), intent(in) :: flann_params
    end function

    !
    ! Adds points to pre-built index.
    !
    ! Params:
    !   index_ptr = pointer to index, must already be built
    !   points = pointer to array of points
    !   rows = number of points to add
    !   columns = feature dimensionality
    !   rebuild_threshold = reallocs index when it grows by factor of
    !     `rebuild_threshold`. A smaller value results is more space efficient
    !     but less computationally efficient. Must be greater than 1.
    !
    ! Returns: 0 if success otherwise -1
    !
    integer(c_int) function flann_add_points(index_id,points,rows,cols,rebuild_threshold) &
                        bind(c,name="flann_add_points")
        import c_int,c_ptr,c_float
        type(c_ptr), value :: index_id
        integer(c_int), intent(in), value :: rows, cols
        real(c_float), intent(in) :: points(cols,rows)
        real(c_float), intent(in), value :: rebuild_threshold
    end function
    integer(c_int) function flann_add_points_float(index_id,points,rows,cols,rebuild_threshold) &
                        bind(c,name="flann_add_points_float")
        import c_int,c_ptr,c_float
        type(c_ptr), value :: index_id
        integer(c_int), intent(in), value :: rows, cols
        real(c_float), intent(in) :: points(cols,rows)
        real(c_float), intent(in), value :: rebuild_threshold
    end function
    integer(c_int) function flann_add_points_double(index_id,points,rows,cols,rebuild_threshold) &
                        bind(c,name="flann_add_points_double")
        import c_int,c_ptr,c_float, c_double
        type(c_ptr), value :: index_id
        integer(c_int), intent(in), value :: rows, cols
        real(c_double), intent(in) :: points(cols,rows)
        real(c_float), intent(in), value :: rebuild_threshold
    end function
    integer(c_int) function flann_add_points_int(index_id,points,rows,cols,rebuild_threshold) &
                        bind(c,name="flann_add_points_int")
        import c_int,c_ptr,c_float
        type(c_ptr), value :: index_id
        integer(c_int), intent(in), value :: rows, cols
        real(c_float), intent(in) :: points(cols,rows)
        real(c_float), intent(in), value :: rebuild_threshold
    end function

    ! 
    ! Removes a point from a pre-built index.
    ! 
    ! index_ptr = pointer to pre-built index.
    ! point_id = index of datapoint to remove.
    !
    integer(c_int) function flann_remove_point(index_ptr,point_id) bind(C,name="flann_remove_point")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
        integer(c_int), intent(in), value :: point_id
    end function
    integer(c_int) function flann_remove_point_float(index_ptr,point_id) bind(C,name="flann_remove_point_float")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
        integer(c_int), intent(in), value :: point_id
    end function
    integer(c_int) function flann_remove_point_double(index_ptr,point_id) bind(C,name="flann_remove_point_double")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
        integer(c_int), intent(in), value :: point_id
    end function
    integer(c_int) function flann_remove_point_int(index_ptr,point_id) bind(C,name="flann_remove_point_int")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
        integer(c_int), intent(in), value :: point_id
    end function

    ! 
    ! Gets a point from a given index position.
    ! 
    ! index_ptr = pointer to pre-built index.
    ! point_id = index of datapoint to get.
    ! 
    ! Returns: pointer to datapoint or NULL on miss
    ! 
    type(c_ptr) function flann_get_point(index_ptr,point_id) bind(C,name="flann_get_point")
        import c_ptr, c_int, c_float
        type(c_ptr), intent(in), value :: index_ptr
        integer(c_int), intent(in), value :: point_id
    end function
    type(c_ptr) function flann_get_point_float(index_ptr,point_id) bind(C,name="flann_get_point_float")
        import c_ptr, c_int, c_float
        type(c_ptr), intent(in), value :: index_ptr
        integer(c_int), intent(in), value :: point_id
    end function
    type(c_ptr) function flann_get_point_double(index_ptr,point_id) bind(C,name="flann_get_point_double")
        import c_ptr, c_int, c_double
        type(c_ptr), intent(in), value :: index_ptr
        integer(c_int), intent(in), value :: point_id
    end function
    type(c_ptr) function flann_get_point_int(index_ptr,point_id) bind(C,name="flann_get_point_int")
        import c_ptr, c_int
        type(c_ptr), intent(in), value :: index_ptr
        integer(c_int), intent(in), value :: point_id
    end function

    ! 
    ! Returns the number of datapoints stored in index.
    ! 
    ! index_ptr = pointer to pre-built index.
    ! 
    integer(c_int) function flann_veclen(index_ptr) bind(C,name="flann_veclen")
        import c_int, c_ptr
        type(c_ptr), intent(in), value :: index_ptr
    end function
    integer(c_int) function flann_veclen_float(index_ptr) bind(C,name="flann_veclen_float")
        import c_int, c_ptr
        type(c_ptr), intent(in), value :: index_ptr
    end function
    integer(c_int) function flann_veclen_double(index_ptr) bind(C,name="flann_veclen_double")
        import c_int, c_ptr
        type(c_ptr), intent(in), value :: index_ptr
    end function
    integer(c_int) function flann_veclen_int(index_ptr) bind(C,name="flann_veclen_int")
        import c_int, c_ptr
        type(c_ptr), intent(in), value :: index_ptr
    end function

    ! 
    ! Returns the dimensionality of datapoints stored in index.
    ! 
    ! index_ptr = pointer to pre-built index.
    !
    integer(c_int) function flann_size(index_ptr) bind(c,name="flann_size")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
    end function
    integer(c_int) function flann_size_float(index_ptr) bind(c,name="flann_size_float")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
    end function
    integer(c_int) function flann_size_double(index_ptr) bind(c,name="flann_size_double")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
    end function
    ! integer(c_int) function flann_size_byte(index_ptr) bind(c,name="flann_size_byte")
    !     import c_int, c_ptr
    !     type(c_ptr), value :: index_ptr
    ! end function
    integer(c_int) function flann_size_int(index_ptr) bind(c,name="flann_size_int")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
    end function

    ! 
    ! Returns the number of bytes consumed by the index.
    ! 
    ! index_ptr = pointer to pre-built index.
    !
    integer(c_int) function flann_used_memory(index_ptr) bind(c,name="flann_used_memory")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
    end function
    integer(c_int) function flann_used_memory_float(index_ptr) bind(c,name="flann_used_memory_float")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
    end function
    integer(c_int) function flann_used_memory_double(index_ptr) bind(c,name="flann_used_memory_double")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
    end function
    ! integer(c_int) function flann_used_memory_byte(index_ptr) bind(c,name="flann_used_memory_byte")
    !     import c_int, c_ptr
    !     type(c_ptr), value :: index_ptr
    ! end function
    integer(c_int) function flann_used_memory_int(index_ptr) bind(c,name="flann_used_memory_int")
        import c_int, c_ptr
        type(c_ptr), value :: index_ptr
    end function

    ! 
    ! Saves the index to a file. Only the index is saved into the file, the dataset corresponding to the index is not saved.
    ! 
    ! @param index_id The index that should be saved
    ! @param filename The filename the index should be saved to
    ! @return Returns 0 on success, negative value on error.
    !  
    integer(c_int) function flann_save_index(index_id,filename) bind(C,name="flann_save_index")
        import c_int,c_ptr,c_char
        type(c_ptr), intent(in), value :: index_id
        character(kind=c_char), intent(in) :: filename(*)
    end function
    integer(c_int) function flann_save_index_float(index_id,filename) bind(C,name="flann_save_index_float")
        import c_int,c_ptr,c_char
        type(c_ptr), intent(in), value :: index_id
        character(kind=c_char), intent(in) :: filename(*)
    end function
    integer(c_int) function flann_save_index_double(index_id,filename) bind(C,name="flann_save_index_double")
        import c_int,c_ptr,c_char
        type(c_ptr), intent(in), value :: index_id
        character(kind=c_char), intent(in) :: filename(*)
    end function
    integer(c_int) function flann_save_index_int(index_id,filename) bind(C,name="flann_save_index_int")
        import c_int,c_ptr,c_char
        type(c_ptr), intent(in), value :: index_id
        character(kind=c_char), intent(in) :: filename(*)
    end function

    ! 
    ! Loads an index from a file.
    ! 
    ! @param filename File to load the index from.
    ! @param dataset The dataset corresponding to the index.
    ! @param rows Dataset tors
    ! @param cols Dataset columns
    ! @return
    ! 
    type(c_ptr) function flann_load_index(filename,dataset,rows,cols) bind(C,name="flann_load_index")
        import c_char,c_int,c_float,c_ptr
        character(kind=c_char), intent(in) :: filename(*)
        integer(c_int), intent(in), value :: rows,cols
        real(c_float), intent(out) :: dataset(cols,rows)
    end function
    type(c_ptr) function flann_load_index_float(filename,dataset,rows,cols) bind(C,name="flann_load_index_float")
        import c_char,c_int,c_float,c_ptr
        character(kind=c_char), intent(in) :: filename(*)
        integer(c_int), intent(in), value :: rows,cols
        real(c_float), intent(out) :: dataset(cols,rows)
    end function
    type(c_ptr) function flann_load_index_double(filename,dataset,rows,cols) bind(C,name="flann_load_index_double")
        import c_char,c_int,c_double,c_ptr
        character(kind=c_char), intent(in) :: filename(*)
        integer(c_int), intent(in), value :: rows,cols
        real(c_double), intent(out) :: dataset(cols,rows)
    end function
    type(c_ptr) function flann_load_index_int(filename,dataset,rows,cols) bind(C,name="flann_load_index_int")
        import c_char,c_int,c_ptr
        character(kind=c_char), intent(in) :: filename(*)
        integer(c_int), intent(in), value :: rows,cols
        integer(c_int), intent(out) :: dataset(cols,rows)
    end function

    ! 
    ! Builds an index and uses it to find nearest neighbors.
    !
    ! Params:
    !  dataset = pointer to a data set stored in row major order
    !  rows = number of rows (features) in the dataset
    !  cols = number of columns in the dataset (feature dimensionality)
    !  testset = pointer to a query set stored in row major order
    !  trows = number of rows (features) in the query dataset (same dimensionality as features in the dataset)
    !  indices = pointer to matrix for the indices of the nearest neighbors of the testset features in the dataset
    !          (must have trows number of rows and nn number of columns)
    !  nn = how many nearest neighbors to return
    !  flann_params = generic flann parameters
    !
    ! Returns: zero or -1 for error
    !
    integer(c_int) function flann_find_nearest_neighbors(dataset,rows,cols,testset,trows,indices,dists,nn,flann_params) &
                    bind(C,name="flann_find_nearest_neighbors")
        import c_int, c_float, FLANNParameters
        integer(c_int), intent(in), value :: rows, cols
        real(c_float), intent(in) :: dataset(cols,rows)
        integer(c_int), intent(in), value :: trows
        real(c_float), intent(in) :: testset(cols,trows)
        integer(c_int), intent(in), value :: nn
        integer(c_int), intent(out) :: indices(nn,trows)
        real(c_float), intent(out) :: dists(nn,trows)
        type(FLANNParameters), intent(in) :: flann_params
    end function
    integer(c_int) function flann_find_nearest_neighbors_float(dataset,rows,cols,testset,trows,indices,dists,nn,flann_params) &
                    bind(C,name="flann_find_nearest_neighbors_float")
        import c_int, c_float, FLANNParameters
        integer(c_int), intent(in), value :: rows, cols
        real(c_float), intent(in) :: dataset(cols,rows)
        integer(c_int), intent(in), value :: trows
        real(c_float), intent(in) :: testset(cols,trows)
        integer(c_int), intent(in), value :: nn
        integer(c_int), intent(out) :: indices(nn,trows)
        real(c_float), intent(out) :: dists(nn,trows)
        type(FLANNParameters), intent(in) :: flann_params
    end function
    integer(c_int) function flann_find_nearest_neighbors_double(dataset,rows,cols,testset,trows,indices,dists,nn,flann_params) &
                    bind(C,name="flann_find_nearest_neighbors_double")
        import c_int, c_double, FLANNParameters
        integer(c_int), intent(in), value :: rows, cols
        real(c_double), intent(in) :: dataset(cols,rows)
        integer(c_int), intent(in), value :: trows
        real(c_double), intent(in) :: testset(cols,trows)
        integer(c_int), intent(in), value :: nn
        integer(c_int), intent(out) :: indices(nn,trows)
        real(c_double), intent(out) :: dists(nn,trows)
        type(FLANNParameters), intent(in) :: flann_params
    end function
    integer(c_int) function flann_find_nearest_neighbors_int(dataset,rows,cols,testset,trows,indices,dists,nn,flann_params) &
                    bind(C,name="flann_find_nearest_neighbors_int")
        import c_int, c_float, FLANNParameters
        integer(c_int), intent(in), value :: rows, cols
        integer(c_int), intent(in) :: dataset(cols,rows)
        integer(c_int), intent(in), value :: trows
        integer(c_int), intent(in) :: testset(cols,trows)
        integer(c_int), intent(in), value :: nn
        integer(c_int), intent(out) :: indices(nn,trows)
        real(c_float), intent(out) :: dists(nn,trows)
        type(FLANNParameters), intent(in) :: flann_params
    end function

    ! 
    ! Searches for nearest neighbors using the index provided
    !
    ! Params:
    !  index_id = the index (constructed previously using flann_build_index).
    !  testset = pointer to a query set stored in row major order
    !  trows = number of rows (features) in the query dataset (same dimensionality as features in the dataset)
    !  indices = pointer to matrix for the indices of the nearest neighbors of the testset features in the dataset
    !          (must have trows number of rows and nn number of columns)
    !  dists = pointer to matrix for the distances of the nearest neighbors of the testset features in the dataset
    !          (must have trows number of rows and 1 column)
    !  nn = how many nearest neighbors to return
    !  flann_params = generic flann parameters
    !
    ! Returns: zero or a number <0 for error
    ! 
    integer(c_int) function flann_find_nearest_neighbors_index(index_id,testset,trows,indices,dists,nn,flann_params) &
                    bind(c,name="flann_find_nearest_neighbors_index")
        import c_int, c_ptr, c_float, FLANNParameters
        type(c_ptr), value :: index_id
        integer(c_int), intent(in), value :: trows
        real(c_float), intent(in) :: testset(*)
        integer(c_int), intent(in), value :: nn
        integer(c_int) :: indices(nn,trows)
        real(c_float) :: dists(nn,trows)  
        type(FLANNParameters) :: flann_params
    end function
    integer(c_int) function flann_find_nearest_neighbors_index_float(index_id,testset,trows,indices,dists,nn,flann_params) &
                    bind(c,name="flann_find_nearest_neighbors_index_float")
        import c_int, c_ptr, c_float, FLANNParameters
        type(c_ptr), value :: index_id
        integer(c_int), intent(in), value :: trows
        real(c_float), intent(in) :: testset(*)
        integer(c_int) :: indices(*)
        real(c_float) :: dists(*)
        integer(c_int), intent(in), value :: nn
        type(FLANNParameters) :: flann_params
    end function
    integer(c_int) function flann_find_nearest_neighbors_index_double(index_id,testset,trows,indices,dists,nn,flann_params) &
                    bind(c,name="flann_find_nearest_neighbors_index_double")
        import c_int, c_ptr, c_float, c_double, FLANNParameters
        type(c_ptr), value :: index_id
        integer(c_int), intent(in), value :: trows
        real(c_double), intent(in) :: testset(*)
        integer(c_int) :: indices(*)
        real(c_double) :: dists(*)
        integer(c_int), intent(in), value :: nn
        type(FLANNParameters) :: flann_params
    end function
    ! integer(c_int) function flann_find_nearest_neighbors_index_byte(index_id,testset,trows,indices,dists,nn,flann_params) &
    !                 bind(c,name="flann_find_nearest_neighbors_index_byte")
    !     import c_int, c_ptr, c_float, FLANNParameters
    !     type(c_ptr), value :: index_id

    !     integer(c_int), intent(in), value :: trows
    !     integer(c_int), intent(out) :: indices(trows,*)
    !     integer(c_int), intent(out) :: dists(trows)
    !     integer(c_int), intent(in), value :: nn
    !     type(FLANNParameters) :: flann_params
    ! end function
    integer(c_int) function flann_find_nearest_neighbors_index_int(index_id,testset,trows,indices,dists,nn,flann_params) &
                    bind(c,name="flann_find_nearest_neighbors_index_int")
        import c_int, c_ptr, c_float, FLANNParameters
        type(c_ptr), value :: index_id
        integer(c_int), intent(in), value :: trows
        integer(c_int), intent(in) :: testset(*)
        integer(c_int) :: indices(*)
        real(c_float) :: dists(*)
        integer(c_int), intent(in), value :: nn
        type(FLANNParameters) :: flann_params
    end function

    ! 
    ! Performs an radius search using an already constructed index.
    ! 
    ! In case of radius search, instead of always returning a predetermined
    ! number of nearest neighbours (for example the 10 nearest neighbours), the
    ! search will return all the neighbours found within a search radius
    ! of the query point.
    ! 
    ! The check parameter in the FLANNParameters below sets the level of approximation
    ! for the search by only visiting "checks" number of features in the index
    ! (the same way as for the KNN search). A lower value for checks will give
    ! a higher search speedup at the cost of potentially not returning all the
    ! neighbours in the specified radius.
    ! 
    ! Returns: count of points within search radius or a number <0 for error
    !
    integer(c_int) function flann_radius_search(index_ptr,query,indices,dists,max_nn,radius,flann_params) &
                    bind(C,name="flann_radius_search")
        import c_int, c_ptr, c_float, FLANNParameters
        type(c_ptr), value :: index_ptr
        real(c_float), intent(in) :: query(*)
        integer(c_int), intent(in), value :: max_nn
        integer(c_int), intent(out) :: indices(max_nn)
        real(c_float), intent(out) :: dists(max_nn)
        real(c_float), intent(in), value :: radius
        type(FLANNParameters), intent(in) :: flann_params
    end function
    integer(c_int) function flann_radius_search_float(index_ptr,query,indices,dists,max_nn,radius,flann_params) &
                    bind(C,name="flann_radius_search_float")
        import c_int, c_ptr, c_float, FLANNParameters
        type(c_ptr), value :: index_ptr
        real(c_float), intent(in) :: query(*)
        integer(c_int), intent(in), value :: max_nn
        integer(c_int), intent(out) :: indices(max_nn)
        real(c_float), intent(out) :: dists(max_nn)
        real(c_float), intent(in), value :: radius
        type(FLANNParameters), intent(in) :: flann_params
    end function
    integer(c_int) function flann_radius_search_double(index_ptr,query,indices,dists,max_nn,radius,flann_params) &
                    bind(C,name="flann_radius_search_double")
        import c_int, c_ptr, c_float, c_double, FLANNParameters
        type(c_ptr), value :: index_ptr
        real(c_double), intent(in) :: query(*)
        integer(c_int), intent(in), value :: max_nn
        integer(c_int), intent(out) :: indices(max_nn)
        real(c_double), intent(out) :: dists(max_nn)
        real(c_float), intent(in), value :: radius
        type(FLANNParameters), intent(in) :: flann_params
    end function
    integer(c_int) function flann_radius_search_int(index_ptr,query,indices,dists,max_nn,radius,flann_params) &
                    bind(C,name="flann_radius_search_int")
        import c_int, c_ptr, c_float, FLANNParameters
        type(c_ptr), value :: index_ptr
        integer(c_int), intent(in) :: query(*)
        integer(c_int), intent(in), value :: max_nn
        integer(c_int), intent(out) :: indices(max_nn)
        integer(c_int), intent(out) :: dists(max_nn)
        real(c_float), intent(in), value :: radius
        type(FLANNParameters), intent(in) :: flann_params
    end function

    ! 
    ! Deletes an index and releases the memory used by it.
    !
    ! Params:
    !  index_id = the index (constructed previously using flann_build_index).
    !  flann_params = generic flann parameters
    !
    ! Returns: zero or a number <0 for error
    !
    integer(c_int) function flann_free_index(index_id,flann_params) bind(c,name="flann_free_index")
        import c_int, c_ptr, FLANNParameters
        type(c_ptr), value :: index_id
        type(FLANNParameters) :: flann_params
    end function
    integer(c_int) function flann_free_index_float(index_id,flann_params) bind(c,name="flann_free_index_float")
        import c_int, c_ptr, FLANNParameters
        type(c_ptr), value :: index_id
        type(FLANNParameters) :: flann_params
    end function
    integer(c_int) function flann_free_index_double(index_id,flann_params) bind(c,name="flann_free_index_double")
        import c_int, c_ptr, FLANNParameters
        type(c_ptr), value :: index_id
        type(FLANNParameters) :: flann_params
    end function
    ! integer(c_int) function flann_free_index_byte(index_id,flann_params) bind(c,name="flann_free_index_byte")
    !     import c_int, c_ptr, FLANNParameters
    !     type(c_ptr), value :: index_id
    !     type(FLANNParameters) :: flann_params
    ! end function
    integer(c_int) function flann_free_index_int(index_id,flann_params) bind(c,name="flann_free_index_int")
        import c_int, c_ptr, FLANNParameters
        type(c_ptr), value :: index_id
        type(FLANNParameters) :: flann_params
    end function

    ! 
    ! Clusters the features in the dataset using a hierarchical kmeans clustering approach.
    ! This is significantly faster than using a flat kmeans clustering for a large number
    ! of clusters.
    !
    ! Params:
    !  dataset = pointer to a data set stored in row major order
    !  rows = number of rows (features) in the dataset
    !  cols = number of columns in the dataset (feature dimensionality)
    !  clusters = number of cluster to compute
    !  result = memory buffer where the output cluster centers are storred
    !  index_params = used to specify the kmeans tree parameters (branching factor, max number of iterations to use)
    !  flann_params = generic flann parameters
    !
    ! Returns: number of clusters computed or a number <0 for error. This number can be different than the number of clusters requested, due to the
    !  way hierarchical clusters are computed. The number of clusters returned will be the highest number of the form
    !  (branch_size-1)*K+1 smaller than the number of clusters requested.
    ! 
    integer(c_int) function flann_compute_cluster_centers(dataset,rows,cols,clusters,result,flann_params) &
                    bind(C,name="flann_compute_cluster_centers")
        import c_int, c_float, FLANNParameters
        integer(c_int), intent(in), value :: rows,cols
        real(c_float), intent(in) :: dataset(cols,rows)
        integer(c_int), intent(in), value :: clusters
        real(c_float), intent(out) :: result(clusters,cols)
        type(FLANNParameters), intent(in) :: flann_params
    end function
    integer(c_int) function flann_compute_cluster_centers_float(dataset,rows,cols,clusters,result,flann_params) &
                    bind(C,name="flann_compute_cluster_centers_float")
        import c_int, c_float, FLANNParameters
        integer(c_int), intent(in), value :: rows,cols
        real(c_float), intent(in) :: dataset(cols,rows)
        integer(c_int), intent(in), value :: clusters
        real(c_float), intent(out) :: result(clusters,cols)
        type(FLANNParameters), intent(in) :: flann_params
    end function
    integer(c_int) function flann_compute_cluster_centers_double(dataset,rows,cols,clusters,result,flann_params) &
                    bind(C,name="flann_compute_cluster_centers_double")
        import c_int, c_double, FLANNParameters
        integer(c_int), intent(in), value :: rows,cols
        real(c_double), intent(in) :: dataset(cols,rows)
        integer(c_int), intent(in), value :: clusters
        real(c_double), intent(out) :: result(clusters,cols)
        type(FLANNParameters), intent(in) :: flann_params
    end function
    integer(c_int) function flann_compute_cluster_centers_int(dataset,rows,cols,clusters,result,flann_params) &
                    bind(C,name="flann_compute_cluster_centers_int")
        import c_int, c_float, FLANNParameters
        integer(c_int), intent(in), value :: rows,cols
        integer(c_int), intent(in) :: dataset(cols,rows)
        integer(c_int), intent(in), value :: clusters
        real(c_float), intent(out) :: result(clusters,cols)
        type(FLANNParameters), intent(in) :: flann_params
    end function

  end interface
    
end module