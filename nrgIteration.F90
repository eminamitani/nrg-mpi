subroutine nrgIteration( iteration )
  use generall
  use parallel
  use hamiltonian, only: numberOfSubspace ! in

  implicit none
  integer, intent(in) :: iteration

  integer :: loadParRank(0:numberOfProcess,1:2)
  logical, allocatable :: ismysubspace(:,:)

  call startCount("nrgIteration")

  if (my_rank .eq. 0) then
     print*, "entering the iteration:", iteration
  end if

  call prepareBasis(iteration)

  allocate( ismysubspace(numberOfSubspace,0:numberOfProcess-1) )
  call loadbalanceDiagonalization2( ismysubspace )

  call diagonalization(iteration,ismysubspace)

  deallocate( ismysubspace )

  !sort the eigenvalue and truncation
  call postprocess(iteration)

  call stopCount
end subroutine nrgIteration
