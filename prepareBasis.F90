subroutine prepareBasis( iteration )
  use generall
  use parallel
  use hamiltonian, only: hami ! inout
  use hamiltonian, only: operation ! in
  use hamiltonian, only: basis_input ! out
  use hamiltonian, only: basis_output ! in
  use hamiltonian, only: numberOfBasis_old ! in
  use hamiltonian, only: numberOfSubspace ! out
  use hamiltonian, only: subspaceInfo ! out
  use hamiltonian, only: allocateType, deallocateType, sortType
  use system, only: allocateType, deallocateType
  implicit none
  integer, intent(in) :: iteration

  integer :: i,j,k,position
  integer :: count,countSpace, eigenvectorIndex
  integer, allocatable ::reference(:)
  character(20) :: format
  character(20) :: Numitr

  call startCount("prepareBasis")

  !definition of the elements in basis_input
  !for basis_input
  !element placed at numberOfVariation+1 : reference : index to define from which eigenstate of last iteration
  !element placed at numberOfVariation+2 : variation : index of eigenstate of the last iteration to distinguish the same conserved quantity
  !element placed at numberOfVariation+3 : operation :what kind of operation is done for creating this basis from the eigenstate of last iteration (add f_{ni\sigma}^dagger etc)

  !at this point, numberOfBasis_old contain the number of eigenstate at last iteration
  hami%numberOfBasis=numberOfBasis_old*hami%numberOfOperation

  if (my_rank .eq. 0) then
     print*, "initial numberOfBasis=", hami%numberOfBasis
  end if

  if( allocated(basis_input) ) then
     call deallocateType(basis_input)
  end if
  call allocateType(basis_input,hami%numberOfBasis,hami%numberOfVariation)
  position=1
  do i=1, numberOfBasis_old
     do j=1, hami%numberOfOperation
        !common operation for type 1 and type 2, just substruct the matrix element.
        !(for make the same behavior in old version of NRG code, here I use substruct)
        !if the Hamiltonian conservatin variable is Parity, not substruct but product is necessary
        basis_input(position)%basis(:) &
             = basis_output(i)%basis(:) - operation(j,:)
        !reference
        basis_input(position)%reference = i
        !variation
        basis_input(position)%variation = basis_output(i)%reference
        !operation
        basis_input(position)%operation = j
        position=position+1
     end do
  end do

  call sortType( basis_input )

  !finish to filling the basis_input and construct subspaceInfo
  !counting the number of subspace
  countSpace=1
  allocate( reference(hami%numberOfVariation) )
  reference(:) = basis_input(1)%basis(:)
  do i=1, hami%numberOfBasis
     if( any(basis_input(i)%basis(:) /= reference(:)) ) then
        reference(:) = basis_input(i)%basis(:)
        countSpace=countSpace+1
     end if
  end do

  deallocate(reference)

  !constructing subspaceInfo
  if (iteration > hami%indexOfFirstIteration) then
     call deallocateType(subspaceInfo)
  end if
  numberOfSubspace=countSpace

  if (my_rank .eq. 0) then
     print*, "numberOfSubspace=", numberOfSubspace
  end if
  call allocateType(subspaceInfo,numberOfSubspace,hami%numberOfVariation)

  count=1
  countSpace=1
  eigenvectorIndex=1

  allocate (reference(hami%numberOfVariation))
  reference(:)=basis_input(1)%basis(:)

  do i=1, hami%numberOfBasis
     if( any(basis_input(i)%basis(:)/=reference(:)) ) then
        subspaceInfo(countSpace)%basis(:)          = reference(:)
        subspaceInfo(countSpace)%dimension         = count-1
        subspaceInfo(countSpace)%count_eigenvector = eigenvectorIndex
        subspaceInfo(countSpace)%start_input       = i-count+1

        eigenvectorIndex=eigenvectorIndex+(count-1)*(count-1)
        count=1
        reference(:) = basis_input(i)%basis(:)
        countSpace=countSpace+1
     end if

     if (i .eq. hami%numberOfBasis) then
        subspaceInfo(countSpace)%basis(:)          = reference(:)
        subspaceInfo(countSpace)%dimension         = count
        subspaceInfo(countSpace)%count_eigenvector = eigenvectorIndex
        subspaceInfo(countSpace)%start_input       = i-count+1
     end if

     count=count+1
  end do

  deallocate(reference)

  call stopCount
end subroutine prepareBasis
