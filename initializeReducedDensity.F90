!****************************************
!Start: 2013/02/28
!E. Minamitani
!this subroutine initialize the reduced density matrix element
!called at the last iteration 1st NRG calculation
!*****************************************

subroutine initializeReducedDensity
  use generall
  use hamiltonian, only: hami ! in
  use hamiltonian, only: numberOfSubspace ! in
  use hamiltonian, only: eigenvalue ! in
  use system, only: fullSubspaceInfo ! inout
  use system, only: pointReducedDensity ! out
  use system, only: reducedDensityMatrix ! out
  use system, only: sizeOfReducedDensity ! in
  use system, only: mapFullSubspaceInfo ! in
  use system, only: fullEigenvalue ! in

  use parallel
  implicit none
  integer :: i,j,k
  integer :: degeneracy
  integer ::indexOfFullSubspaceInfo
  integer ::numberOfBasisInSub
  integer :: isubl, isubr
  double precision :: partitionFunction, eigenval

  call startCount("iniRedDen")
  !simple 0K case
  !count the degeneracy of the ground state
  degeneracy=0
  do i=1, hami%numberOfBasis
     if(eigenvalue(i) < thresouldOfDegeneracy) then
        degeneracy=degeneracy+1
     end if
  end do
  partitionFunction=dble(degeneracy)

  if (my_rank .eq. 0) then
     print*, "initialize reduced density"
     print*, "reduced density matrix size:", sizeOfReducedDensity
     print*, "last iteration :: partition function=", partitionFunction
  end if

  allocate( reducedDensityMatrix(1:sizeOfReducedDensity) )
  reducedDensityMatrix=0.0d0
  indexOfFullSubspaceInfo=mapFullSubspaceInfo(hami%numberOfIteration)
  pointReducedDensity=1

  !fill up for last iteraion part
  do i=1, numberOfSubspace
     !fullSubspaceInfo
     !numberOfVariation+1 : the degeneracy of this subspace before truncation
     !numberOfVariation+2 : degeneracy of this subspace after truncation
     !numberOfVariation+3 : the starting point of this subspace in fullEigenvector
     !numberOfVariation+4 : the starting point of this subspace in fullBasisInput
     !numberOfVariation+5 : the starting point of this subspace in fullBasisOutput (& fullEigenvalue)
     !numberOfVariation+6 : the starting point of this subspace in reducedDensityMatrix, filled in the calculation of reduced density matrix part
     fullSubspaceInfo(indexOfFullSubspaceInfo+i-1)%start_matrix = pointReducedDensity
     numberOfBasisInSub=fullSubspaceInfo(indexOfFullSubspaceInfo+i-1)%dimension_after

     do isubl=1, numberOfBasisInSub
        do isubr=1, numberOfBasisInSub
           if(isubl .eq. isubr) then
              eigenval=fullEigenvalue(fullSubspaceInfo(indexOfFullSubspaceInfo+i-1)%start_output+isubl-1)
              if (eigenval < thresouldOfDegeneracy) then
                 reducedDensityMatrix(pointReducedDensity) = 1.0d0/partitionFunction
              end if

              pointReducedDensity=pointReducedDensity+1
           else
              reducedDensityMatrix(pointReducedDensity) =0.0d0
              pointReducedDensity=pointReducedDensity+1
           end if
        end do
     end do
  end do

  call stopCount
end subroutine initializeReducedDensity
