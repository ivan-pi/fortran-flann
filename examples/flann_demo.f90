program flann_demo

  use iso_c_binding
  use flann
  
  implicit none

  integer, parameter :: nn = 3
  integer, parameter :: npoints = 9000
  integer, parameter :: dim = 128

  integer, parameter :: qcount = 1000

  real(c_float), allocatable :: dataset(:,:)
  real(c_float), allocatable :: query(:,:)

  type(flann_index) :: idx

  integer, allocatable :: indices(:,:)
  real(c_float), allocatable :: dists(:,:)

  call flann_log_verbosity(FLANN_LOG_DEBUG)

  allocate(dataset(dim,npoints))
  allocate(query(dim,qcount))

  call random_number(dataset)
  call random_number(query)

  !! Construct a randomized kd-tree index using 4 kd-trees
  idx = flann_index(dataset,KDTreeSingleIndexParams(),'Euclidean')

  call idx%build_index()
  
  allocate(indices(nn,qcount))
  allocate(dists(nn,qcount))

  indices = 0
  dists = 0

  !! Perform a knn searh, using 128 checks
  call idx%knn_search(query,indices,dists,nn,search_params(checks=128))

end program
