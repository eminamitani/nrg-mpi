module generall
  use count_module

  double precision ,parameter :: pi=3.141592653589793238462643383279502884197d0
  double precision ,parameter :: numericalThresould=1.0D-12
  double precision, parameter :: eigenvectorThresould =1.0D-14
  double precision, parameter :: thresouldOfDegeneracy=1.0D-8

  !support interface to reallocate the 2-dimensional array
  interface reallocate
     module procedure reallocate_1D_real
  end interface


contains
  subroutine sortWithIndex( size, array, indeces )
    integer, intent(in) :: size
    double precision, intent(inout) :: array(size)
    integer, intent(out) :: indeces(size)

    integer :: i, j
    double precision :: swap_value
    integer :: swap_index
    
    do i=1, size
       indeces(i) = i
    end do

    do i=1+1, size
       j=i
       do while(j>1)
          if( array(j) < array(j-1) ) then
             swap_value = array(j)
             array(j) = array(j-1)
             array(j-1) = swap_value

             swap_index = indeces(j)
             indeces(j) = indeces(j-1)
             indeces(j-1) = swap_index
          end if
          j=j-1
       end do
    end do
  end subroutine sortWithIndex

  subroutine sortByIndex( size, array, indeces )
    integer, intent(in) :: size
    double precision, intent(inout) :: array(size)
    integer, intent(inout) :: indeces(size)

    integer :: i, j
    double precision :: swap_value
    integer :: swap_index
    
    do i=1+1, size
       j=i
       do while(j>1)
          if( indeces(j) < indeces(j-1) ) then
             swap_value = array(j)
             array(j) = array(j-1)
             array(j-1) = swap_value

             swap_index = indeces(j)
             indeces(j) = indeces(j-1)
             indeces(j-1) = swap_index
          end if
          j=j-1
       end do
    end do
  end subroutine sortByIndex

  subroutine reallocate_1D_real(array1D, resize)
    implicit none
    double precision, allocatable::array1D(:)
    integer :: arraySize
    integer, intent(in) :: resize
    double precision, allocatable :: tmp(:)

    arraySize=size(array1D)

    allocate(tmp(arraySize))

    tmp=array1D

    deallocate(array1D)
    allocate(array1D(1:resize))
    array1D=0.0d0
    array1D(1:arraySize)=tmp(1:arraySize)
    deallocate(tmp)
  end subroutine reallocate_1D_real
end module generall

module parallel
  integer :: numberOfProcess
  integer :: my_rank
end module parallel


module io_param
  character(20) :: willson_chain_diagonal="diagonal"
  character(20) :: willson_chain_nondiagonal="nondiagonal"
  integer, parameter :: errorlog=1001
  character(20) :: errorfile="error.log"
  integer,parameter :: buglog=2001
  character(50) :: bugfile
  integer,parameter :: eigenodd=2002
  character(20) :: fileEigenOdd="eigenodd.txt"
  integer,parameter :: eigeneven=2003
  character(20) :: fileEigenEven="eigeneven.txt"
  integer, parameter :: stateInfomation=2004
  character(20) :: fileInfo = "keepingState.txt"

  integer, parameter :: thermoXn=3001
  integer,parameter ::thermoSn=3002
  integer, parameter ::thermoCn=3003
  character(20)::fileXn="Xn"
  character(20)::fileSn="Sn"
  character(20)::fileCn="Cn"
end module io_param

!!-------------------------------------------------------------------
module SubspaceInfo_module
  implicit none
  type SubspaceInfo_type
     integer, pointer :: memory(:)   => null() ! memory for array
     integer, pointer :: basis(:)    => null() ! 1:numberOfVariation
     integer, pointer :: dimension   => null() ! numberOfVariation+1
     integer, pointer :: count_eigenvector => null() ! numberOfVariation+2
     integer, pointer :: start_input => null() ! numberOfVariation+3
  end type SubspaceInfo_type

  interface allocateType
     module procedure allocateSubspaceInfoArray
  end interface
  interface deallocateType
     module procedure deallocateSubspaceInfoArray
  end interface
  interface reallocateType
     module procedure reallocateSubspaceInfoArray
  end interface
  interface sizeofType
     module procedure sizeofSubspaceInfoArray
  end interface
  interface locationType
     module procedure locationSubspaceInfoArray
  end interface
  interface copyType
     module procedure copySubspaceInfo
     module procedure copySubspaceInfoArray
  end interface

contains

  integer function sizeofSubspaceInfoArray( array ) result(s)
    type(SubspaceInfo_type), intent(in) :: array(:)
    integer :: n
    s = size(array)*(size(array(1)%basis) + 3)*sizeof(n)
  end function sizeofSubspaceInfoArray

  function locationSubspaceInfoArray( array ) result(location)
    type(SubspaceInfo_type), intent(in) :: array(:)
    integer, pointer :: location

    location => array(1)%basis(1)
  end function locationSubspaceInfoArray
  
  subroutine allocateSubspaceInfoArray( array, num, len )
    type(SubspaceInfo_type), allocatable :: array(:)
    integer, intent(in) :: num, len

    integer :: i, k
    integer, pointer :: memory(:)

    if( allocated(array) ) then
       call deallocateSubspaceInfoArray( array )
    end if

    allocate( array(num) )
    allocate( memory(num*(len+3)) )
 
    array(1)%memory => memory
    k=1
    do i=1, num
       array(i)%basis => memory(k:k+len-1)
       k=k+len
       array(i)%dimension => memory(k)
       k=k+1
       array(i)%count_eigenvector => memory(k)
       k=k+1
       array(i)%start_input => memory(k)
       k=k+1
    end do

  end subroutine allocateSubspaceInfoArray

  subroutine deallocateSubspaceInfoArray( array )
    type(SubspaceInfo_type), allocatable :: array(:)    
    integer :: i

    if( .not. allocated(array) ) return

    deallocate(array(1)%memory)
    do i=1, size(array)
       array(i)%basis => null()
       array(i)%dimension => null()
       array(i)%count_eigenvector => null()
       array(i)%start_input => null()
    end do

    deallocate(array)
  end subroutine deallocateSubspaceInfoArray

  subroutine reallocateSubspaceInfoArray( array, resize, len )
    type(SubspaceInfo_type), allocatable :: array(:)
    integer, intent(in) :: resize, len

    type(SubspaceInfo_type), allocatable :: work(:)

    call allocateSubspaceInfoArray( work, size(array), len )
    call copySubspaceInfoArray( work, array )
    call deallocateSubspaceInfoArray( array )
    call allocateSubspaceInfoArray( array, resize, len )
    call copySubspaceInfoArray( array, work )
    call deallocateSubspaceInfoArray( work )
  end subroutine reallocateSubspaceInfoArray


  subroutine copySubspaceInfo( variableB, variableA )
    type(SubspaceInfo_type), intent(inout) :: variableB
    type(SubspaceInfo_type), intent(in)    :: variableA

    variableB%basis(:)          = variableA%basis(:)
    variableB%dimension         = variableA%dimension
    variableB%count_eigenvector = variableA%count_eigenvector
    variableB%start_input       = variableA%start_input
  end subroutine copySubspaceInfo

  subroutine copySubspaceInfoArray( arrayB, arrayA )
    type(SubspaceInfo_type), intent(inout) :: arrayB(:)
    type(SubspaceInfo_type), intent(in)    :: arrayA(:)

    integer :: i
    do i=1, size(arrayA)
       call copySubspaceInfo( arrayB(i), arrayA(i) )
    end do
  end subroutine copySubspaceInfoArray

end module SubspaceInfo_module

           !fullSubspaceInfo
           !numberOfVariation+1 : the degeneracy of this subspace before truncation
           !numberOfVariation+2 : degeneracy of this subspace after truncation
           !numberOfVariation+3 : the starting point of this subspace in fullEigenvector
           !numberOfVariation+4 : the starting point of this subspace in fullBasisInput
           !numberOfVariation+5 : the starting point of this subspace in fullBasisOutput (& fullEigenvalue)
           !numberOfVariation+6 : the starting point of this subspace in reducedDensityMatrix, filled in the calculation of reduced density matrix part

!!-------------------------------------------------------------------
module FullSubspaceInfo_module
  implicit none
  type FullSubspaceInfo_type
     integer, pointer :: memory(:)    => null() ! memory for array
     integer, pointer :: basis(:)     => null() ! 1:numberOfVariation
     integer, pointer :: dimension_before => null() ! numberOfVariation+1
     integer, pointer :: dimension_after  => null() ! numberOfVariation+2
     integer, pointer :: start_eigen  => null() ! numberOfVariation+3
     integer, pointer :: start_input  => null() ! numberOfVariation+4
     integer, pointer :: start_output => null() ! numberOfVariation+5
     integer, pointer :: start_matrix => null() ! numberOfVariation+6
  end type FullSubspaceInfo_type

  interface allocateType
     module procedure allocateFullSubspaceInfoArray
  end interface
  interface deallocateType
     module procedure deallocateFullSubspaceInfoArray
  end interface
  interface reallocateType
     module procedure reallocateFullSubspaceInfoArray
  end interface
  interface sizeofType
     module procedure sizeofFullSubspaceInfoArray
  end interface
  interface locationType
     module procedure locationFullSubspaceInfoArray
  end interface
  interface copyType
     module procedure copyFullSubspaceInfo
     module procedure copyFullSubspaceInfoArray
  end interface

contains

  subroutine allocateFullSubspaceInfoArray( array, num, len )
    type(FullSubspaceInfo_type), allocatable :: array(:)
    integer, intent(in) :: num, len

    integer :: i, k
    integer, pointer :: memory(:)

    if( allocated(array) ) then
       call deallocateFullSubspaceInfoArray( array )
    end if

    allocate( array(num) )
    allocate( memory(num*(len+6)) )

    array(1)%memory => memory
    k=1
    do i=1, num
       array(i)%basis => memory(k:k+len-1)
       k=k+len
       array(i)%dimension_before => memory(k)
       k=k+1
       array(i)%dimension_after => memory(k)
       k=k+1
       array(i)%start_eigen  => memory(k)
       k=k+1
       array(i)%start_input  => memory(k)
       k=k+1
       array(i)%start_output => memory(k)
       k=k+1
       array(i)%start_matrix => memory(k)
       k=k+1
    end do

  end subroutine allocateFullSubspaceInfoArray

  subroutine deallocateFullSubspaceInfoArray( array )
    type(FullSubspaceInfo_type), allocatable :: array(:)    
    integer :: i

    if( .not. allocated(array) ) return

    deallocate(array(1)%memory)
    do i=1, size(array)
       array(i)%basis => null()
       array(i)%dimension_before => null()
       array(i)%dimension_after =>  null()
       array(i)%start_eigen  =>  null()
       array(i)%start_input  =>  null()
       array(i)%start_output =>  null()
       array(i)%start_matrix =>  null()
    end do

    deallocate(array)
  end subroutine deallocateFullSubspaceInfoArray

  subroutine reallocateFullSubspaceInfoArray( array, resize, len )
    type(FullSubspaceInfo_type), allocatable :: array(:)
    integer, intent(in) :: resize, len

    type(FullSubspaceInfo_type), allocatable :: work(:)

    call allocateFullSubspaceInfoArray( work, size(array), len )
    call copyFullSubspaceInfoArray( work, array )
    call deallocateFullSubspaceInfoArray( array )
    call allocateFullSubspaceInfoArray( array, resize, len )
    call copyFullSubspaceInfoArray( array, work )
    call deallocateFullSubspaceInfoArray( work )
  end subroutine reallocateFullSubspaceInfoArray

  integer function sizeofFullSubspaceInfoArray( array ) result(s)
    type(FullSubspaceInfo_type), intent(in) :: array(:)
    integer :: n
    s = size(array)*(size(array(1)%basis) + 6)*sizeof(n)
  end function sizeofFullSubspaceInfoArray

  function locationFullSubspaceInfoArray( array ) result(location)
    type(FullSubspaceInfo_type), intent(in) :: array(:)
    integer, pointer :: location

    location => array(1)%basis(1)
  end function locationFullSubspaceInfoArray

  subroutine copyFullSubspaceInfo( variableB, variableA )
    type(FullSubspaceInfo_type), intent(inout) :: variableB
    type(FullSubspaceInfo_type), intent(in)  :: variableA

    variableB%basis(:)          = variableA%basis(:)
    variableB%dimension_before  = variableA%dimension_before
    variableB%dimension_after   = variableA%dimension_after
    variableB%start_eigen       = variableA%start_eigen
    variableB%start_input       = variableA%start_input
    variableB%start_output      = variableA%start_output
    variableB%start_matrix      = variableA%start_matrix
  end subroutine copyFullSubspaceInfo

  subroutine copyFullSubspaceInfoArray( arrayB, arrayA )
    type(FullSubspaceInfo_type), intent(inout) :: arrayB(:)
    type(FullSubspaceInfo_type), intent(in)    :: arrayA(:)

    integer :: i
    do i=1, size(arrayA)
       call copyFullSubspaceInfo( arrayB(i), arrayA(i) )
    end do
  end subroutine copyFullSubspaceInfoArray

end module FullSubspaceInfo_module



!!-------------------------------------------------------------------
module Basis_module
  use count_module
  implicit none

  type Basis_type
     integer, pointer :: memory(:) => null() ! memory for array
     integer, pointer :: basis(:)  => null() ! 1:numberOfVariation
     integer, pointer :: reference => null() ! numberOfVariation+1
     integer, pointer :: variation => null() ! numberOfVariation+2
     integer, pointer :: operation => null() ! numberOfVariation+3
  end type Basis_type

  interface allocateType
     module procedure allocateBasisArray
  end interface
  interface deallocateType
     module procedure deallocateBasisArray
  end interface
  interface reallocateType
     module procedure reallocateBasisArray
  end interface
  interface sizeofType
     module procedure sizeofBasisArray
  end interface
  interface locationType
     module procedure locationBasisArray
  end interface
  interface copyType
     module procedure copyBasis
     module procedure copyBasisArray
  end interface
  interface sortType
     module procedure sortBasisArray
  end interface

contains
  subroutine allocateBasisArray( array, num, len )
    type(Basis_type), allocatable :: array(:)
    integer, intent(in) :: num, len

    integer :: i, k
    integer, pointer :: memory(:)

    if( allocated(array) ) then
       call deallocateBasisArray( array )
    end if

    allocate( array(num) )
    allocate( memory(num*(len+3)) )

    array(1)%memory => memory
    k=1
    do i=1, num
       array(i)%basis => memory(k:k+len-1)
       k=k+len
       array(i)%reference => memory(k)
       k=k+1
       array(i)%variation => memory(k)
       k=k+1
       array(i)%operation => memory(k)
       k=k+1
    end do
  end subroutine allocateBasisArray

  subroutine deallocateBasisArray( array )
    type(Basis_type), allocatable :: array(:)    
    integer :: i

    if( .not. allocated(array) ) return

    deallocate(array(1)%memory)
    do i=1, size(array)
       array(i)%basis     => null()
       array(i)%reference => null()
       array(i)%variation => null()
       array(i)%operation => null()
    end do
    deallocate(array)
  end subroutine deallocateBasisArray

  subroutine reallocateBasisArray( array, resize, len )
    type(Basis_type), allocatable :: array(:)
    integer, intent(in) :: resize, len

    type(Basis_type), allocatable :: work(:)

    call allocateBasisArray( work, size(array), len )
    call copyBasisArray( work, array )
    call deallocateBasisArray( array )
    call allocateBasisArray( array, resize, len )
    call copyBasisArray( array, work )
    call deallocateBasisArray( work )
  end subroutine reallocateBasisArray

  integer function sizeofBasisArray( array ) result(s)
    type(Basis_type), intent(in) :: array(:)
    integer :: n
    s = size(array)*(size(array(1)%basis) + 3)*sizeof(n)
  end function sizeofBasisArray

  function locationBasisArray( array ) result(location)
    type(Basis_type), intent(in) :: array(:)
    integer, pointer :: location

    location => array(1)%basis(1)
  end function locationBasisArray

  subroutine copyBasis( variableB, variableA )
    type(Basis_type), intent(inout) :: variableB
    type(Basis_type), intent(in)    :: variableA

    variableB%basis(:)  = variableA%basis(:)
    variableB%reference = variableA%reference
    variableB%variation = variableA%variation
    variableB%operation = variableA%operation
  end subroutine copyBasis

  subroutine copyBasisArray( arrayB, arrayA )
    type(Basis_type), intent(inout) :: arrayB(:)
    type(Basis_type), intent(in)    :: arrayA(:)

    integer :: i
    do i=1, size(arrayA)
       call copyBasis( arrayB(i), arrayA(i) )
    end do
  end subroutine copyBasisArray

  subroutine sortBasisArray( array )
    type(Basis_type), intent(inout) :: array(:)

    integer :: n, m

    call startCount("sortBasis")
    do n=1, size(array)
       do m=n+1, size(array)
          if( isgreaterBasis( array(n), array(m) ) ) then
             call swapBasis( array(n), array(m) )
          end if
       end do
    end do

    call stopCount
    return

  contains
    logical function isgreaterBasis( variableA, variableB ) result(isgreater)
      type(Basis_type), intent(inout) :: variableA
      type(Basis_type), intent(inout) :: variableB

      integer :: i
      do i=1, size(variableA%basis)
         if( variableA%basis(i) > variableB%basis(i) ) then
            isgreater = .true.
            return 
         else if( variableA%basis(i) < variableB%basis(i) ) then
            isgreater = .false.
            return 
         else !! variableA%basis(i)==variableB%basis(i)
            ! try next icol
         end if
      end do

      if( variableA%reference > variableB%reference ) then
         isgreater = .true.
         return 
      else if( variableA%reference < variableB%reference ) then
         isgreater = .false.
         return 
      end if

      isgreater = .false.
      !! variableA%basis(:)==variableB%basis(:) .and.
      !! variableA%reference==variableB%reference
    end function isgreaterBasis

    subroutine swapBasis( variableA, variableB )
      type(Basis_type), intent(inout) :: variableA
      type(Basis_type), intent(inout) :: variableB

      integer :: i

      do i=1, size(variableA%basis)
         call swapInteger( variableA%basis(i), variableB%basis(i) )
      end do
      call swapInteger( variableA%reference, variableB%reference )
      call swapInteger( variableA%variation, variableB%variation )
      call swapInteger( variableA%operation, variableB%operation )
    end subroutine swapBasis

    subroutine swapInteger( a, b )
      integer, intent(inout) :: a, b
      integer :: w

      w = a
      a = b
      b = w
    end subroutine swapInteger

  end subroutine sortBasisArray

end module Basis_module


module system
  use FullSubspaceInfo_module
  use Basis_module
  use generall

  type Parameter_type
     double precision :: bandWidth
     double precision :: ztrick
     double precision :: Lambda
     ! define the scaling method
     ! 1. Wilson original type, 2. Campo's method
     integer :: scalingType
     ! reference number of trancation of state
     integer :: truncation
     ! determine whether the cfs type of calculation is done or not
     logical :: cfsnrg
     logical :: noimpurity
     logical :: loadRestart
     logical :: saveRestart
  end type Parameter_type
  type(Parameter_type) :: param

  double precision, parameter :: truncationThreshold=0.000001d0
  double precision, parameter :: truncationLimitScale=1.2d0

  !integer :: cfsNRGstep
  integer, allocatable :: mapFullSubspaceInfo(:)
  integer :: pointFullBasisInput, pointFullBasisOutput
  integer :: pointFullSubspace, pointReducedDensity
  !change to 64 bit integer to parse large size matrix
  integer(kind=8) :: pointFullEigenvector
  integer :: elementNumFullBasisInput, elementNumFullBasisOutput
  integer :: elementNumFullSubspace, elementNumReducedDensity, elementNumFullEigenVector

  type(FullSubspaceInfo_type), allocatable :: fullSubspaceInfo(:)
  type(Basis_type), allocatable :: fullBasisOutput(:)
  type(Basis_type), allocatable :: fullBasisInput(:)
  double precision, allocatable :: fullEigenvector(:)
  double precision, allocatable :: fullEigenvalue(:)

  integer :: sizeOfReducedDensity
  double precision, allocatable :: reducedDensityMatrix(:)
  integer :: columnOfFullSubspaceInfo

  !the scale factor to calculate the thermodynamical quantities
  double precision :: betabar=0.8


contains
  subroutine read_parameter
    double precision :: bandWidth
    double precision :: ztrick
    double precision :: Lambda
    integer :: scalingType
    integer :: truncation
    logical :: cfsnrg
    logical :: noimpurity
    logical :: loadRestart, saveRestart
    namelist /nrg/ bandWidth, ztrick, Lambda, scalingType, &
         truncation, cfsnrg, noimpurity, loadRestart, saveRestart

    ! default values
    cfsnrg      = .false.
    noimpurity  = .false.
    loadRestart = .false.
    saveRestart = .false.

    open(10,file="parameter.txt")
    read(10,nrg)
    close(10)

    param%bandWidth   = bandWidth
    param%ztrick      = ztrick
    param%Lambda      = Lambda
    param%scalingType = scalingType
    param%truncation  = truncation
    param%cfsnrg      = cfsnrg
    param%noimpurity  = noimpurity
    param%loadRestart = loadRestart
    param%saveRestart = saveRestart

  end subroutine read_parameter

  double precision function scale(N)
    implicit none
    integer :: N

    select case (param%scalingType)
    case(0)
       !no scaling
       scale=param%Lambda**(-dble(N)*0.5d0)

    case(1)
       !return the energy scale of each iteration
       !in original Wilson method
       scale=param%bandWidth*0.5d0*(1.0d0+1.0d0/param%Lambda) &
            *param%Lambda**(-dble(N)*0.5d0-param%ztrick+1.0d0)

    case(2)
       !return the energy scale of each iteration
       !in Campo+z-trick( see NRG Lijunbrura)
       scale=param%Lambda**(-dble(N)*0.5d0-param%ztrick) &
            *param%bandWidth*(param%Lambda-1.0d0)/log(param%Lambda)
    end select
  end function scale
end module system

module hamiltonian
  use generall
  use parallel
  use SubspaceInfo_module
  use Basis_module
  !this module keeps model hamiltonian specific set up
  implicit none

  type HamiltonianInfo_type
     ! number of Wilson chain: channel number
     integer :: numberOfIteration
     ! number of conserved physical quantities
     integer :: numberOfChain
     ! contain the starting and ending element of charge related
     ! coservation value and spin Sz related conservation value
     integer :: numberOfVariation
     ! number of type of Invariant matrix element for localized
     ! electron excitation spectrum
     integer :: conservationBlock(2,2)
     ! number of type of invariant matrix for conduction electron part
     integer :: numberOfMatrix
     ! contain number of basis in each iteration
     integer :: numberOfConductionMatrix
     ! number of coefficient type for invariant matrix element
     ! for conduction electron and impurity site each.
     integer :: numberOfBasis
     ! number of coefficient type for invariant matrix element
     ! for conduction electron and impurity site each.
     integer :: numberOfCoefficient
     ! number of variation of the electronic state adding to
     ! each iteration
     integer :: numberOfCoefficientImpurity
     integer :: numberOfOperation
     !number of BCS paring coefficient
     integer :: numberOfBCSCoefficient
     integer :: indexOfFirstIteration
     ! whether calculate the spectrum related value
     logical :: flag_spectrum
     ! if other local operator
     logical :: flag_local
     !these entities are for CFS-NRG run
     integer :: startTrancation
     logical :: flag_truncation
     !whether BCS type superconducting lead exist or not
     logical :: flag_BCS

  end type HamiltonianInfo_type
  type(HamiltonianInfo_type) :: hami

  !integer :: typeOfHamiltonian
  !array for contain the hybridization matrix element between impurity and conduction electron
  !usual case, I set numberOfHybridization=numberOfMatrix (spin dependent hybridization is also possible)
  !for define the size of this array, I set one integer named numberOfHybridization
  integer :: numberOfHybridization
  !array for contain the hybridization between impurity and conduction chain. dimension size =number of impurity*2 * number of chain
  ! double precision, allocatable :: hybridization(:,:)
  !array for wilson chain. diagonal part and nondiagonal part
  double precision, allocatable :: chain_diagonal(:,:)
  double precision, allocatable :: chain_nondiagonal(:,:)
  !additional chain element in the case of BCS superconductor
  double precision,allocatable :: chain_BCS(:,:)

  integer ::lowOfInput, lowOfOutput, lowOfSubspace
  !contain number of basis in each iteration
  integer :: numberOfBasis_old,numberOfBasis_full


  !the law value of ground state in each iteration
  double precision :: groundStateEnergy
  type(SubspaceInfo_type), allocatable :: subspaceInfo(:)
  type(Basis_type), allocatable :: basis_input(:)
  type(Basis_type), allocatable :: basis_output(:)

  integer, allocatable :: operation(:,:)
  integer, allocatable :: coefficient_invariant_matrix_type(:,:)
  double precision,allocatable :: coefficient_invariant_matrix(:)
  !order of elements in coefficient_invariant_matrix_type :: operation_right, operation_left, typeofmatrix
  !file "coefficient.dat"
  !operation_right, operation_left, typeofmatrix, coefficient

  double precision, allocatable :: coefficient_invariant_matrix_spectrum(:,:)
  !the coefficient required to calculate the single particle exitation spectrum
  !since in the calculation of the matrix element above -1 iteration,
  !only the matrix element between the basis derived from same operation
  !have nonzero value.
  !Thus, the coefficient_invariant_matrix _spectrum becomes one dimensional array
  !and the index is corresponds to the kind of operation.
  integer,allocatable :: conservation_difference_spectrum(:,:)
  !in the calculation of the invatriant matrix for spectrum, use this array to find the pair of the subspace which can possess matrix elements

  integer,allocatable :: coefficient_diagonalization_type(:,:)
  double precision, allocatable :: coefficient_diagonalization(:)
  !define the basis pair contribute to BCS type paring
  integer, allocatable :: coefficient_BCS(:,:)

  double precision, allocatable :: invariant_matrix(:,:,:)
  double precision, allocatable :: invariant_matrix_spectrum(:,:,:)
  double precision, allocatable :: invariant_matrix_spectrum_initial(:,:,:)
  integer :: numberOfSubspace

  double precision, allocatable :: eigenvector(:)
  integer :: sizeOfEigenvector

  double precision, allocatable :: eigenvalue(:)

  !for set up MPI communication of type:eigenvalue
  integer :: iblock(2), idisp(2), itype(2)

  !for keeping the information at initial state
  integer :: numberOfSubspace_initial, sizeOfEigenvector_initial, numberOfBasis_initial

contains
  subroutine setupHamiltonian
    include 'mpif.h'
    integer ::ierr

    integer :: numberOfIteration
    integer :: numberOfChain
    integer :: numberOfVariation
    integer :: conservationBlock(2,2)
    integer :: numberOfMatrix
    integer :: numberOfConductionMatrix
    integer :: numberOfBasis
    integer :: numberOfCoefficient
    integer :: numberOfCoefficientImpurity
    integer :: numberOfOperation
    integer :: numberOfBCSCoefficient
    integer :: indexOfFirstIteration
    logical :: flag_spectrum
    logical :: flag_local
    logical :: flag_BCS

    namelist /hamiltonianInfo/ &
         numberOfIteration, numberOfChain, numberOfVariation, conservationBlock, &
         numberOfMatrix, numberOfConductionMatrix, numberOfBasis, numberOfCoefficient, numberOfBCSCoefficient, &
         numberOfCoefficientImpurity, numberOfOperation, indexOfFirstIteration, &
         flag_spectrum, flag_local, flag_BCS

    call startCount("setupHami")

    if (my_rank .eq. 0) then
       ! default values
       numberOfCoefficientImpurity = 0
       flag_BCS=.false.
       numberOfBCSCoefficient=0

       open(10,file="hamiltonianInfo.dat")
       read(10,hamiltonianInfo)
       close(10)

       hami%numberOfIteration           = numberOfIteration
       hami%numberOfChain               = numberOfChain
       hami%numberOfVariation           = numberOfVariation
       hami%conservationBlock(:,:)      = conservationBlock(:,:)
       hami%numberOfMatrix              = numberOfMatrix
       hami%numberOfConductionMatrix    = numberOfConductionMatrix
       hami%numberOfBasis               = numberOfBasis
       hami%numberOfCoefficient         = numberOfCoefficient
       hami%numberOfCoefficientImpurity = numberOfCoefficientImpurity
       hami%numberOfOperation           = numberOfOperation
       hami%numberOfBCSCoefficient      = numberOfBCSCoefficient
       hami%indexOfFirstIteration       = indexOfFirstIteration
       hami%flag_spectrum               = flag_spectrum
       hami%flag_local                  = flag_local
       hami%flag_BCS                    = flag_BCS
       hami%startTrancation             = hami%numberOfIteration+1
       hami%flag_truncation             = .false.
    end if

    call MPI_BCAST(hami,sizeof(hami), MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
!!$    call MPI_BCAST(numberOfIteration, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!!$    call MPI_BCAST(numberOfChain,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$    call MPI_BCAST(numberOfVariation,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$    call MPI_BCAST(conservationBlock(1,1),size(conservationBlock),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$    call MPI_BCAST(numberOfMatrix,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$    call MPI_BCAST(numberOfConductionMatrix,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$    call MPI_BCAST(numberOfBasis,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$    call MPI_BCAST(numberOfCoefficient,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$    call MPI_BCAST(numberOfCoefficientImpurity,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$    call MPI_BCAST(numberOfOperation,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$    call MPI_BCAST(indexOfFirstIteration, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
!!$    call MPI_BCAST(flag_spectrum,1,MPI_LOGICAL,0,MPI_COMM_WORLD, ierr)
!!$    call MPI_BCAST(flag_local,1,MPI_LOGICAL,0,MPI_COMM_WORLD, ierr)

    lowOfInput=hami%numberOfVariation+3
    lowOfOutput=hami%numberOfVariation+1
    lowOfSubspace=hami%numberOfVariation+3

    !SubspaceInfo explanation
    !numberOfVariation+1 :: number of the basis belongs to the subspace, corresponds to rmax in old version code
    !numberOfVariation+2 :: start position of the eigenvector corresponds to the subspace
    !numberOfVariation+3 :: start position of BasisInput corresponds to the subspace

    call allocateType( basis_output, &
         hami%numberOfBasis, hami%numberOfVariation )
    allocate( eigenvalue(hami%numberOfBasis) )
    allocate( operation(hami%numberOfOperation, &
         hami%numberOfVariation+hami%numberOfChain) )
    allocate( coefficient_invariant_matrix_type &
         ( hami%numberOfCoefficient,3) )
    allocate( coefficient_invariant_matrix(hami%numberOfCoefficient) )
    allocate( coefficient_diagonalization_type &
         ( hami%numberOfCoefficient,4) )
    allocate( coefficient_diagonalization(hami%numberOfCoefficient) )
    allocate( invariant_matrix( hami%numberOfBasis, &
         hami%numberOfBasis, hami%numberOfConductionMatrix) )

    if(hami%flag_spectrum) then
       allocate( coefficient_invariant_matrix_spectrum &
            ( hami%numberOfMatrix, hami%numberOfOperation) )
       !initialize invariant_matrix_spectrum
       if (my_rank .eq. 0) then
          print*, "invariant_matrix_spectrum size=", &
               hami%numberOfBasis*hami%numberOfBasis*hami%numberOfMatrix
       end if
       allocate( invariant_matrix_spectrum &
            ( hami%numberOfBasis, hami%numberOfBasis, hami%numberOfMatrix) )
       allocate( invariant_matrix_spectrum_initial &
            ( hami%numberOfBasis, hami%numberOfBasis, hami%numberOfMatrix) )
       allocate( conservation_difference_spectrum &
            ( hami%numberOfMatrix, hami%numberOfVariation) )
    end if

    if (hami%flag_BCS) then
        allocate(coefficient_BCS(1:hami%numberOfBCSCoefficient,1:3))
    end if

    !allocate(difference(1:numberOfDifference, 1:numberOfVariation))

    call stopCount
  end subroutine setupHamiltonian

end module hamiltonian




module spectrum
  use generall
  use hamiltonian, only: Basis_type ! type
  use hamiltonian, only: SubspaceInfo_type ! type
  implicit none

  !container for spectrum (reterded Green's function calculation) matrix element.
  !since this constructs discrete spectrum often called as spike.
  !the number of spike itself is huge and unpredictable, thus, I binned it after calculation
  !pair of energy and value, positive energy--> excitation negativeenergy--> deexcitation
  double precision, allocatable :: binnedSpikePositive(:,:),binnedSpikeNegative(:,:)

  integer :: numberOfBin !the number of scales in binnedSpike, the number of element if numberOfBin*numberOfMatrix
  integer, parameter ::meshOfBin=1000

  !the number of element of binnedSpike is determined the energy scale in NRG run
  double precision :: Emax, log10Emax
  !following the 2mipurity former code
  double precision, parameter:: Emin=1.0D-8
  double precision, parameter:: log10Emin=-15
  integer, allocatable :: positiveSpectrumOperation(:,:), negativeSpectrumOperation(:,:)

  !to reduce calculation, restore the results of 1st NRG run
  double precision, allocatable :: restoreEigenvalue(:)
  double precision, allocatable :: restoreEigenVector(:)
  double precision, allocatable :: restore_invariant_matrix_spectrum(:)
  type(Basis_type), allocatable :: restore_basis_output(:)

  type(Basis_type), allocatable :: restore_basis_input(:)
  type(SubspaceInfo_type), allocatable :: restore_subspaceInfo(:)

  !for restore the datas in 1st NRG run for CFS spectrum calculation
  type restore
     integer :: eigenvalue_start
     integer :: eigenvalue_number
     integer :: eigenvector_start
     integer :: eigenvector_number
     integer :: subspace_start
     integer :: subspace_number
     integer :: basis_output_start
     integer :: invariant_matrix_spectrum_start
     integer :: invariant_matrix_spectrum_number
     integer :: numberOfBasis
  end type restore

  type(restore), allocatable :: restoreInfo(:)

contains

  subroutine setSpectrumInfo
    use hamiltonian, only: hami ! in
!!$    use spectrum, only: positiveSpectrumOperation ! out
!!$    use spectrum, only: negativeSpectrumOperation ! out
    call startCount("setSpect")

    !TODO fill
    allocate( positiveSpectrumOperation(hami%numberOfMatrix,hami%numberOfVariation) )
    allocate( negativeSpectrumOperation(hami%numberOfMatrix,hami%numberOfVariation) )

    call stopCount
  end subroutine setSpectrumInfo

! CFS calculation
  subroutine spectrumElementDGEMM &
       ( dim_left,rmax_left, eigenvec_min_left, basis_min_left, basis_max_left, &
         dim_right,rmax_right, eigenvec_min_right, basis_min_right, basis_max_right, &
         matrixType, spectrumArray )
    use hamiltonian, only: basis_input ! in
    use hamiltonian, only: eigenvector ! in
    use hamiltonian, only: invariant_matrix_spectrum ! in
    use hamiltonian, only: coefficient_invariant_matrix_spectrum ! in

    integer, intent(in) :: dim_left, dim_right
    integer, intent(in) :: rmax_left, rmax_right
    integer, intent(in) :: basis_min_left, basis_min_right
    integer, intent(in) :: basis_max_left, basis_max_right
    integer, intent(in) :: eigenvec_min_left, eigenvec_min_right
    integer, intent(in) :: matrixType
    double precision, intent(out) :: spectrumArray(dim_left,dim_right)

    integer :: endLeft, endRight
    integer :: operation_left, operation_right
    integer :: reference_left, reference_right
    integer :: startLeft2, startRight2

    integer :: il, ir, isubl, isubr
    double precision, allocatable :: vectorL(:,:)
    double precision, allocatable :: vectorR(:,:)
    double precision, allocatable :: coefMatrix(:,:)
    double precision, allocatable :: tmp1(:,:)
    double precision :: spike_temp, coef
    external :: dgemm
    double precision, parameter :: alpha=1.0d0, beta=0.0d0
    integer:: rowA, columnA, rowB, columnB, rowC, columnC

    call startCount("spectElem")


    allocate( vectorL (dim_left, rmax_left ) )
    allocate( vectorR (rmax_right, dim_right) )

    vectorL=0.0d0
    vectorR=0.0d0
    spectrumArray=0.0d0

    allocate (tmp1(dim_left,rmax_right ))
    tmp1=0.0d0

    allocate (coefMatrix(rmax_left,rmax_right))
    coefMatrix=0.0d0

    do isubl=1, dim_left
       vectorL(isubl, 1:rmax_left)=eigenvector(eigenvec_min_left+(isubl-1)*rmax_left:eigenvec_min_left+isubl+rmax_left-1)
    end do

    do isubr=1, dim_right
       vectorR(1:rmax_right, isubr)=eigenvector(eigenvec_min_right+(isubr-1)*rmax_right:eigenvec_min_right+isubr+rmax_right-1)
    end do

         do isubl=basis_min_left, basis_max_left

            operation_left=basis_input(isubl)%operation
            reference_left=basis_input(isubl)%reference

            loop_sub_r: do isubr=basis_min_right, basis_max_right

                operation_right=basis_input(isubr)%operation
                reference_right=basis_input(isubr)%reference

                if (operation_right .ne. operation_left) then
                    cycle loop_sub_r
                end if



                    coefMatrix(isubl-basis_min_left+1,isubr-basis_min_right+1) = coefficient_invariant_matrix_spectrum &
                        (matrixType,operation_left)*invariant_matrix_spectrum(reference_right,reference_left,matrixType)


            end do loop_sub_r
        end do

            rowA=dim_left
            columnA=rmax_left
            rowB=rmax_left
            columnB=rmax_right
            rowC=dim_left
            columnC=rmax_right

            call dgemm('N','N',rowA,columnB,columnA,alpha,vectorL,rowA,&
                coefMatrix,rowB,beta,tmp1,rowC)

            rowA=dim_left
            columnA=rmax_right
            !vectorR-->B
            rowB=rmax_right
            columnB=dim_right
            !invariant -->C
            rowC=dim_left
            columnC=dim_right

            !tmp1*vectorR^T
            call dgemm('N','N',rowA, columnB,columnA,alpha,tmp1,rowA,&
                vectorR,rowB,beta,spectrumArray,rowC)

    deallocate( vectorL  )
    deallocate( vectorR )
    deallocate (tmp1)
    deallocate (coefMatrix)

    call stopCount



  end subroutine spectrumElementDGEMM

  double precision function spectrumElement &
       ( startLeft, startRight, &
       maxVariationLeft, maxVariationRight, &
       variationLeft, variationRight, &
       startEigenvecLeft, startEigenvecRight, matrixType )
    use hamiltonian, only: basis_input ! in
    use hamiltonian, only: eigenvector ! in
    use hamiltonian, only: invariant_matrix_spectrum ! in
    use hamiltonian, only: coefficient_invariant_matrix_spectrum ! in

    integer, intent(in) :: startLeft, startRight
    integer, intent(in) :: maxVariationLeft, maxVariationRight
    integer, intent(in) :: variationLeft, variationRight
    integer, intent(in) :: startEigenvecLeft, startEigenvecRight
    integer, intent(in) :: matrixType

    integer :: endLeft, endRight
    integer :: operationLeft, operationRight
    integer :: referenceLeft, referenceRight
    integer :: startLeft2, startRight2

    integer :: il, ir
    double precision, allocatable :: eigenvector_left(:)
    double precision, allocatable :: eigenvector_right(:)
    double precision :: spike_temp, coef

    call startCount("spectElem")

    endLeft  = startLeft  + maxVariationLeft  - 1
    endRight = startRight + maxVariationRight - 1
    startLeft2  = startEigenvecLeft &
         - startLeft + maxVariationLeft*(variationLeft-1)
    startRight2 = startEigenvecRight &
         - startRight + maxVariationRight*(variationRight-1)

    allocate( eigenvector_left (startLeft :endLeft ) )
    allocate( eigenvector_right(startRight:endRight) )

    eigenvector_left(startLeft:endLeft) &
         = eigenvector(startLeft2+startLeft:startLeft2+endLeft)
    eigenvector_right(startRight:endRight) &
         = eigenvector(startRight2+startRight:startRight2+endRight)

    spectrumElement = 0.0d0
    do il=startLeft, endLeft !minimum index of basis_input to maximum index of basis_input
       if( eigenvector_left(il) == 0.0d0 ) cycle
       operationLeft = basis_input(il)%operation
       referenceLeft = basis_input(il)%reference !reference the index (correspond to the order of Out element)

       spike_temp = 0.0d0
       do ir=startRight, endRight
          if( eigenvector_right(ir) == 0.0d0 ) cycle
          operationRight = basis_input(ir)%operation
          if( operationLeft .ne. operationRight ) cycle
          referenceRight = basis_input(ir)%reference

!!$          spectrumElement = spectrumElement &
!!$               + coefficient_invariant_matrix_spectrum &
!!$               ( matrixType, operationLeft ) &
!!$               * eigenvector_left(il) &
!!$               * eigenvector_right(ir) &
!!$               * invariant_matrix_spectrum &
!!$               (referenceRight,referenceLeft,matrixType)
          spike_temp = spike_temp &
               + eigenvector_right(ir) &
               * invariant_matrix_spectrum &
               (referenceRight,referenceLeft,matrixType)
       end do ! ir
       if( spike_temp == 0.0d0 ) cycle

       spectrumElement = spectrumElement &
            + coefficient_invariant_matrix_spectrum &
            ( matrixType, operationLeft ) &
            * eigenvector_left(il) * spike_temp
    end do ! il

    deallocate( eigenvector_left  )
    deallocate( eigenvector_right )

    call stopCount
  end function spectrumElement

  !subroutine for binning the spike
  subroutine binningSpike(container, spike, matrixType)
    integer, intent(in) :: matrixType
    !container is binnedSpikeNegative pr binnedSpikePositive
    double precision, allocatable :: container(:,:)
    !single spike to bin
    double precision :: spike(1,1:2)
    integer :: i1, i2, positionReference
    double precision :: logreference,weightBin, position

    call startCount("binSpike")

    i1=0
    i2=0
    logreference=log10(abs(spike(1,1)))
    position=(logreference-log10Emin)*meshOfBin+1.0d0
    !using "int" in fortran makes truncate under dicimal point
    positionReference =int(position)

    if (position  .le. 0.0d0) then
       i1=1
       i2=1
    else if (position .ge. numberOfBin) then
       i1=numberOfBin
       i2=numberOfBin
    end if

    if (i1==i2 .and. i1>0 .and. i2>0) then
       container(i1+(matrixType-1)*numberOfBin,2) &
            = container(i1+(matrixType-1)*numberOfBin,2)+spike(1,2)
    else
       i1=positionReference
       i2=positionReference+1
       weightBin=position-dble(positionReference)
       container(i1+(matrixType-1)*numberOfBin,2) &
            = container(i1+(matrixType-1)*numberOfBin,2)+spike(1,2)*(1.0d0-weightBin)
       container(i2+(matrixType-1)*numberOfBin,2) &
            = container(i2+(matrixType-1)*numberOfBin,2)+spike(1,2)*weightBin
    end if

    call stopCount
  end subroutine binningSpike

  subroutine findSubspace(tmpVariation,isub,flag)
    use hamiltonian, only: hami ! in
    use hamiltonian, only: subspaceInfo ! in
    use hamiltonian, only: numberOfSubspace ! in

    integer, allocatable :: tmpVariation(:)
    logical :: flag
    integer :: isub
    integer :: i, iv, diff

    call startCount("findSub")

    flag=.false.

    do i=1, numberOfSubspace
       diff=0
       do iv=1, hami%numberOfVariation
          diff=diff+abs(subspaceInfo(i)%basis(iv)-tmpVariation(iv))
       end do

       if(diff .eq. 0) then
          isub=i
          flag=.true.
          exit
       end if
    end do

    call stopCount
  end subroutine findSubspace

  subroutine findFullSubspace(tmpVariation,isub, itr, flag)
    use hamiltonian, only: hami ! in
    use system, only: fullSubspaceInfo ! in
    use system, only: mapFullSubspaceInfo ! in

    integer, allocatable :: tmpVariation(:)
    integer :: isub
    integer:: i,iv,diff, itr
    logical :: flag

    call startCount("findFull")

    flag=.false.

    do i=mapFullSubspaceInfo(itr), mapFullSubspaceInfo(itr+1)-1
       diff=0
       do iv=1, hami%numberOfVariation
          diff=diff+abs(fullSubspaceInfo(i)%basis(iv)-tmpVariation(iv))
       end do

       if(diff .eq. 0) then
          flag=.true.
          isub=i
          exit
       end if
    end do

    call stopCount
  end subroutine findFullSubspace

  integer function findStartOutput(tmpVariation)
    use hamiltonian, only: hami ! in
    use hamiltonian, only: basis_output ! in

    integer, allocatable :: tmpVariation(:)
    integer :: i, iv, diff

    call startCount("findStart")

    do i=1, hami%numberOfBasis
       diff=0
       do iv=1, hami%numberOfVariation
          diff=diff+abs(basis_output(i)%basis(iv)-tmpVariation(iv))
       end do

       if(diff .eq. 0) then
          findStartOutput=i
          exit
       end if
    end do

    call stopCount
  end function findStartOutput



  !for reduce the calculational cost in spectrum calculation, store the results in 1st NEG calculation
  !and reuse them in the second calculation

  subroutine initRestore
    use generall
    use parallel
    use hamiltonian, only: hami ! in
    use hamiltonian, only: numberOfBasis_full ! in
    use hamiltonian, only: numberOfSubspace ! in
    use hamiltonian, only: eigenvector ! in
    use hamiltonian, only: eigenvalue ! in
    use hamiltonian, only: basis_input ! in
    use hamiltonian, only: basis_output ! in
    use hamiltonian, only: invariant_matrix_spectrum ! in
    use hamiltonian, only: subspaceInfo ! in
    use hamiltonian, only: allocateType, copyType ! routine
    use system, only: allocateType, copyType ! routine
!!$    use spectrum, only: restoreInfo ! out
!!$    use spectrum, only: restore_subspaceInfo ! out
!!$    use spectrum, only: restore_basis_output ! out
!!$    use spectrum, only: restore_basis_input ! out
!!$    use spectrum, only: restoreEigenvalue ! out
!!$    use spectrum, only: restoreEigenvector ! out
!!$    use spectrum, only: restore_invariant_matrix_spectrum ! out
    implicit none
    integer :: i,j, imat, icount
    integer ::eigenvectorSize
    !this subrouitne is to initialize the container of restore information in 1st NRG calculation
    !this should be called in the postprocess of 1st iteration

    call startCount("iniRestore")

    eigenvectorSize=size(eigenvector)
    allocate( restoreInfo(hami%indexOfFirstIteration:hami%numberOfIteration+1) )
    allocate( restoreEigenvalue(numberOfBasis_full) )
    !tentative allocation
    allocate( restoreEigenVector(eigenvectorSize) )
    allocate( restore_invariant_matrix_spectrum &
         ( hami%numberOfBasis*hami%numberOfBasis*hami%numberOfMatrix) )
    call allocateType( restore_subspaceInfo, &
         numberOfSubspace, hami%numberOfVariation )
    call allocateType( restore_basis_output, &
         numberOfBasis_full, hami%numberOfVariation )
    call allocateType( restore_basis_input, &
         numberOfBasis_full, hami%numberOfVariation )

    restoreInfo(hami%indexOfFirstIteration)%eigenvalue_start = 1
    restoreInfo(hami%indexOfFirstIteration)%eigenvalue_number = numberOfBasis_full
    restoreInfo(hami%indexOfFirstIteration)%basis_output_start = 1
    restoreInfo(hami%indexOfFirstIteration)%eigenvector_start = 1
    restoreInfo(hami%indexOfFirstIteration)%eigenvector_number = eigenvectorSize
    restoreInfo(hami%indexOfFirstIteration)%subspace_start = 1
    restoreInfo(hami%indexOfFirstIteration)%subspace_number = numberOfSubspace
    restoreInfo(hami%indexOfFirstIteration)%invariant_matrix_spectrum_start = 1
    restoreInfo(hami%indexOfFirstIteration)%invariant_matrix_spectrum_number = &
         hami%numberOfMatrix*hami%numberOfBasis*hami%numberOfBasis
    restoreInfo(hami%indexOfFirstIteration)%numberOfBasis = hami%numberOfBasis

    do i=1, numberOfBasis_full
       call copyType( restore_basis_input(i), basis_input(i) )
       call copyType( restore_basis_output(i), basis_output(i) )
       restoreEigenvalue(i) = eigenvalue(i)
    end do

    icount=1
    do imat=1, hami%numberOfMatrix
       do j=1, hami%numberOfBasis
          do i=1, hami%numberOfBasis
             restore_invariant_matrix_spectrum(icount)=invariant_matrix_spectrum(i,j,imat)
             icount=icount+1
          end do
       end do
    end do

    do i=1, eigenvectorSize
       restoreEigenVector(i)=eigenvector(i)
    end do

    do i=1, numberOfSubspace
       call copyType( restore_subspaceInfo(i), subspaceInfo(i) )
    end do

    if (my_rank .eq. 0) then
       print*, "indexOfFirstIteration=", hami%indexOfFirstIteration
       print*, "initial information of restore, numberOfBasis_full, numberOfBasis=", numberOfBasis_full, hami%numberOfBasis
       print*, "restoreInfo%eigenvalue_number", &
            restoreInfo(hami%indexOfFirstIteration)%eigenvalue_number
    end if

    call stopCount
  end subroutine initRestore


  subroutine appendRestore( iteration )
    use generall
    use parallel
    use hamiltonian, only: hami ! in
    use hamiltonian, only: numberOfBasis_full ! in
    use hamiltonian, only: numberOfSubspace ! in
    use hamiltonian, only: eigenvector ! in
    use hamiltonian, only: eigenvalue ! in
    use hamiltonian, only: basis_input ! in
    use hamiltonian, only: basis_output ! in
    use hamiltonian, only: invariant_matrix_spectrum ! in
    use hamiltonian, only: subspaceInfo ! in
    use hamiltonian, only: reallocateType, copyType ! routine
    use system, only: reallocateType, copyType ! routine
!!$    use spectrum, only: restoreInfo ! inout
!!$    use spectrum, only: restoreEigenvalue ! inout
!!$    use spectrum, only: restoreEigenvector ! inout
!!$    use spectrum, only: restore_invariant_matrix_spectrum ! inout
!!$    use spectrum, only: restore_basis_output ! inout
!!$    use spectrum, only: restore_basis_input ! inout
!!$    use spectrum, only: restore_subspaceInfo ! inout
    implicit none
    integer, intent(in) :: iteration

    integer :: i, j, imat, icount
    integer :: eigenvalue_start_last, eigenvalue_number_last
    integer :: eigenvalue_start_present, eigenvalue_end
    integer :: eigenvectorSize, eigenvector_start_last
    integer :: eigenvector_number_last, eigenvector_start_present, eigenvector_end
    integer :: subspace_start_last, subspace_number_last, subspace_start_present, subspace_end
    integer :: basis_output_start_last, basis_output_start_present, basis_output_end
    integer :: invariant_matrix_spectrum_start_last, invariant_matrix_spectrum_start_present
    integer :: invariant_matrix_spectrum_end
    integer :: invariant_matrix_spectrum_number_last

    call startCount("appRestore")

    !restore eigenvalue
    eigenvalue_start_last = restoreInfo(iteration-1)%eigenvalue_start
    eigenvalue_number_last = restoreInfo(iteration-1)%eigenvalue_number
    eigenvalue_start_present = eigenvalue_start_last+eigenvalue_number_last
    restoreInfo(iteration)%eigenvalue_start = eigenvalue_start_present
    restoreInfo(iteration)%eigenvalue_number = numberOfBasis_full
    eigenvalue_end = eigenvalue_start_present+numberOfBasis_full-1
    call reallocate( restoreEigenvalue, eigenvalue_end )

    if (my_rank .eq. 0) then
       print*, "this iteration is ", iteration
       print*, "eigenvalue_number_last=",eigenvalue_number_last
    end if

    !restore basis_output
    basis_output_start_last = restoreInfo(iteration-1)%basis_output_start
    basis_output_start_present = basis_output_start_last+eigenvalue_number_last
    basis_output_end = basis_output_start_present+numberOfBasis_full-1
    restoreInfo(iteration)%basis_output_start = basis_output_start_present
    call reallocateType( restore_basis_output, basis_output_end, hami%numberOfVariation )
    call reallocateType( restore_basis_input, basis_output_end, hami%numberOfVariation )

    do i=1, numberOfBasis_full
       restoreEigenvalue(eigenvalue_start_present+i-1) = eigenvalue(i)
       call copyType( restore_basis_output(basis_output_start_present+i-1), basis_output(i) )
       call copyType( restore_basis_input(basis_output_start_present+i-1), basis_input(i) )
    end do

    !restore eigenvector
    eigenvectorSize = size(eigenvector)
    eigenvector_start_last = restoreInfo(iteration-1)%eigenvector_start
    eigenvector_number_last = restoreInfo(iteration-1)%eigenvector_number
    eigenvector_start_present = eigenvector_start_last+eigenvector_number_last
    eigenvector_end = eigenvector_start_present+eigenvectorSize-1
    restoreInfo(iteration)%eigenvector_start = eigenvector_start_present
    restoreInfo(iteration)%eigenvector_number = eigenvectorSize
    call reallocate( restoreEigenVector, eigenvector_end )

    do i=1, eigenvectorSize
       restoreEigenVector(eigenvector_start_present+i-1) = eigenvector(i)
    end do

    !restore subspaceInfo
    subspace_start_last = restoreInfo(iteration-1)%subspace_start
    subspace_number_last = restoreInfo(iteration-1)%subspace_number
    subspace_start_present = subspace_start_last+subspace_number_last
    subspace_end = subspace_start_present + numberOfSubspace - 1
    restoreInfo(iteration)%subspace_start = subspace_start_present
    restoreInfo(iteration)%subspace_number = numberOfSubspace

    if (my_rank .eq. 0) then
       print*, "numberOfSubspace,", numberOfSubspace
       print*, "subspace_end,", subspace_end
    end if

    call reallocateType( restore_subspaceInfo, &
         subspace_end, hami%numberOfVariation )
    do i=1, numberOfSubspace
       call copyType &
            ( restore_subspaceInfo(subspace_start_present+i-1), &
            subspaceInfo(i) )
    end do

    !restore invariant Matrix spectrum
    if(iteration<hami%numberOfIteration) then
       invariant_matrix_spectrum_start_last &
            = restoreInfo(iteration-1)%invariant_matrix_spectrum_start
       invariant_matrix_spectrum_number_last &
            = restoreInfo(iteration-1)%invariant_matrix_spectrum_number
       invariant_matrix_spectrum_start_present &
            = invariant_matrix_spectrum_start_last &
            + invariant_matrix_spectrum_number_last
       invariant_matrix_spectrum_end &
            = invariant_matrix_spectrum_start_present &
            +hami%numberOfBasis*hami%numberOfBasis*hami%numberOfMatrix-1
       restoreInfo(iteration)%invariant_matrix_spectrum_start &
            = invariant_matrix_spectrum_start_present
       restoreInfo(iteration)%invariant_matrix_spectrum_number &
            = hami%numberOfBasis*hami%numberOfBasis*hami%numberOfMatrix
       call reallocate(restore_invariant_matrix_spectrum,invariant_matrix_spectrum_end)

       icount=1
       do imat=1, hami%numberOfMatrix
          do j=1, hami%numberOfBasis
             do i=1, hami%numberOfBasis
                restore_invariant_matrix_spectrum &
                     (invariant_matrix_spectrum_start_present+icount-1) &
                     = invariant_matrix_spectrum(i,j,imat)
                icount=icount+1
             end do
          end do
       end do

    end if

    restoreInfo(iteration)%numberOfBasis=hami%numberOfBasis

    call stopCount
  end subroutine appendRestore


  subroutine getRestore( iteration )
    use generall
    use parallel
    use hamiltonian, only: hami ! inout
    use hamiltonian, only: numberOfBasis_full ! out
    use hamiltonian, only: numberOfBasis_old ! out
    use hamiltonian, only: numberOfSubspace ! out
    use hamiltonian, only: eigenvector ! out
    use hamiltonian, only: eigenvalue ! out
    use hamiltonian, only: basis_input ! out
    use hamiltonian, only: basis_output ! out
    use hamiltonian, only: invariant_matrix_spectrum ! out
    use hamiltonian, only: invariant_matrix_spectrum_initial ! in
    use hamiltonian, only: subspaceInfo ! out
    use hamiltonian, only: allocateType, deallocateType, copyType ! routine
    use system, only: allocateType, deallocateType, copyType ! routine
!!$    use spectrum, only: restoreInfo ! in
!!$    use spectrum, only: restoreEigenvalue ! in
!!$    use spectrum, only: restoreEigenvector ! in
!!$    use spectrum, only: restore_invariant_matrix_spectrum ! in
!!$    use spectrum, only: restore_basis_output ! in
!!$    use spectrum, only: restore_basis_input ! in
!!$    use spectrum, only: restore_subspaceInfo ! in
    implicit none
    integer, intent(in) :: iteration

    integer :: eigenvalue_start, eigenvalue_number, eigenvector_start, eigenvector_number
    integer :: basis_output_start, invariant_matrix_spectrum_start,invariant_matrix_spectrum_number
    integer :: subspace_start, subspace_number
    integer ::i, j, imat, icount
    character(20) :: Numitr
    character(50) :: format

    call startCount("getRestore")

    eigenvalue_start = restoreInfo(iteration)%eigenvalue_start
    eigenvalue_number = restoreInfo(iteration)%eigenvalue_number
    eigenvector_start = restoreInfo(iteration)%eigenvector_start
    eigenvector_number = restoreInfo(iteration)%eigenvector_number
    basis_output_start = restoreInfo(iteration)%basis_output_start

    subspace_start = restoreInfo(iteration)%subspace_start
    subspace_number = restoreInfo(iteration)%subspace_number

    numberOfSubspace = subspace_number
    hami%numberOfBasis = restoreInfo(iteration)%numberOfBasis
    numberOfBasis_full = restoreInfo(iteration)%eigenvalue_number
    numberOfBasis_old = restoreInfo(iteration-1)%numberOfBasis
    if (my_rank .eq. 0) then
       print*, "numberOfBasis, numberOfBasis_full, numberOfBasis_old=", hami%numberOfBasis, numberOfBasis_full, numberOfBasis_old
    end if

    if( allocated(eigenvalue) ) deallocate(eigenvalue)
    if( allocated(eigenvector) ) deallocate(eigenvector)
    call deallocateType(basis_output)
    if( allocated(invariant_matrix_spectrum) ) deallocate(invariant_matrix_spectrum)
    call deallocateType(subspaceInfo)
    call deallocateType(basis_input)

    if (iteration > hami%indexOfFirstIteration) then
       invariant_matrix_spectrum_start = restoreInfo(iteration-1)%invariant_matrix_spectrum_start
       invariant_matrix_spectrum_number = restoreInfo(iteration-1)%invariant_matrix_spectrum_number
    else
       invariant_matrix_spectrum_number = size(invariant_matrix_spectrum_initial)
    end if

    allocate( eigenvalue(eigenvalue_number) )
    allocate( eigenvector(eigenvector_number) )
    call allocateType( basis_output, eigenvalue_number, hami%numberOfVariation )
    call allocateType( subspaceInfo, subspace_number, hami%numberOfVariation )

    allocate( invariant_matrix_spectrum &
         ( numberOfBasis_old, numberOfBasis_old, hami%numberOfMatrix) )

    call allocateType( basis_input, eigenvalue_number, hami%numberOfVariation )

    do i=1, eigenvalue_number
       eigenvalue(i) = restoreEigenvalue(eigenvalue_start+i-1)
       call copyType( basis_output(i), restore_basis_output(basis_output_start+i-1) )
       call copyType( basis_input(i), restore_basis_input(basis_output_start+i-1) )
    end do

    do i=1, eigenvector_number
       eigenvector(i) = restoreEigenVector(eigenvector_start+i-1)
    end do

    do i=1, subspace_number
       call copyType( subspaceInfo(i), &
            restore_subspaceInfo(subspace_start+i-1) )
    end do

    if(iteration>hami%indexOfFirstIteration ) then
       icount=1
       do imat=1, hami%numberOfMatrix
          do j=1, numberOfBasis_old
             do i=1, numberOfBasis_old
                invariant_matrix_spectrum(i,j,imat) &
                     = restore_invariant_matrix_spectrum(invariant_matrix_spectrum_start+icount-1)
                icount=icount+1
             end do
          end do
       end do
    else !in the case, the truncation is started at initial iteration
       invariant_matrix_spectrum(:,:,:) &
            = invariant_matrix_spectrum_initial(:,:,:)
    end if

    call stopCount
  end subroutine getRestore

  subroutine cfsSpectrumCalculation( iteration )
    use generall
    use parallel
    use hamiltonian, only: hami ! in
    use hamiltonian, only: numberOfBasis_full ! in
    implicit none

    integer, intent(in) :: iteration

    call startCount("cfsSpectrum")

    if (my_rank .eq. 0) then
       print*, "iteration=", iteration
    end if
    call getRestore( iteration )

    if (iteration .eq. hami%numberOfIteration) then
       call spikeZeroTemperature( iteration )
    else
       call spikeCFS( iteration, numberOfBasis_full)
    end if

    call stopCount
  end subroutine cfsSpectrumCalculation

end module spectrum

module restart
  use generall
  use parallel
  use system
  use hamiltonian
  use spectrum
  implicit none

  integer, parameter :: MODE_RESTART_NON = 1000
  integer, parameter :: MODE_RESTART_NRG = 1001
  integer, parameter :: MODE_RESTART_CFS = 1002

contains
  subroutine saveRestart( iteration )
    include 'mpif.h'
    integer :: ierr

    integer, intent(in) :: iteration
    integer, parameter :: fd = 100

    if( my_rank == 0 ) then
       write(*,*) "saving the latest CFS data to restart.dat"
!!$       open(fd,file="restart.dat",status='replace',form='binary')
       open(fd,file="restart.dat",status='replace',form="unformatted",access="stream")
       !!-------
       write(fd) MODE_RESTART_CFS
       write(fd) iteration
       !!-------
       write(fd) hami
       !!-------
       write(fd) lbound(restoreInfo)
       write(fd) ubound(restoreInfo)
       write(fd) restoreInfo(:)
       !!-------
       write(fd) size(restore_subspaceInfo)
       write(fd) restore_subspaceInfo(1)%memory(:)
       !!-------
       write(fd) size(restore_basis_input)
       write(fd) restore_basis_input(1)%memory(:)
       !!-------
       write(fd) size(restore_basis_output)
       write(fd) restore_basis_output(1)%memory(:)
       !!-------
       write(fd) size(restore_invariant_matrix_spectrum)
       write(fd) restore_invariant_matrix_spectrum(:)
       !!-------
       write(fd) size(invariant_matrix_spectrum_initial,1)
       write(fd) size(invariant_matrix_spectrum_initial,2)
       write(fd) size(invariant_matrix_spectrum_initial,3)
       write(fd) invariant_matrix_spectrum_initial(:,:,:)
       !!-------
       write(fd) size(restoreEigenvalue)
       write(fd) restoreEigenvalue(:)
       !!-------
       write(fd) size(restoreEigenVector)
       write(fd) restoreEigenVector(:)
       !!-------
       write(fd) size(reducedDensityMatrix)
       write(fd) reducedDensityMatrix(:)
       !!-------
       write(fd) size(fullSubspaceInfo)
       write(fd) fullSubspaceInfo(1)%memory(:)
       !!-------
       write(fd) lbound(mapFullSubspaceInfo)
       write(fd) ubound(mapFullSubspaceInfo)
       write(fd) mapFullSubspaceInfo(:)
       !!-------
       write(fd) size(binnedSpikePositive)
       write(fd) binnedSpikePositive(:,:)
       !!-------
       write(fd) size(binnedSpikeNegative)
       write(fd) binnedSpikeNegative(:,:)
       !!-------
       close(fd)
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine saveRestart

  subroutine loadRestart( iteration_restart, mode_restart )
    include 'mpif.h'
    integer :: ierr

    integer, intent(out) :: iteration_restart
    integer, intent(out) :: mode_restart

    logical :: ex
    integer, parameter :: fd = 100
    integer :: arraysize
    integer :: arraysize1,arraysize2,arraysize3

    inquire(file="restart.dat",exist=ex)
    if( .not. ex ) then
      if( my_rank == 0 ) then
        print*, "restart file not found"
      end if
      call MPI_Abort( MPI_COMM_WORLD, 0, ierr )
    end if

    write(*,*) "loading the last CFS data from restart.dat"
!!$    open(fd,file="restart.dat",status='old',form='binary')
    open(fd,file="restart.dat",status='old',form="unformatted",access="stream")
    !!-------
    read(fd) mode_restart
    if( mode_restart /= MODE_RESTART_CFS ) then
      if( my_rank == 0 ) then
        print*, "broken restart file"
      end if
      call MPI_Abort( MPI_COMM_WORLD, 0, ierr )
    end if
    !!-------
    read(fd) iteration_restart
    !!-------
    read(fd) hami
    !!-------
    read(fd) arraysize1,arraysize2
    if( allocated(restoreInfo) ) deallocate(restoreInfo)
    allocate( restoreInfo(arraysize1:arraysize2) )
    read(fd) restoreInfo(:)
    !!-------
    read(fd) arraysize
    call allocateType( restore_subspaceInfo,arraysize,hami%numberOfVariation )
    read(fd) restore_subspaceInfo(1)%memory(:)
    !!-------
    read(fd) arraysize
    call allocateType( restore_basis_input,arraysize,hami%numberOfVariation )
    read(fd) restore_basis_input(1)%memory(:)
    !!-------
    read(fd) arraysize
    call allocateType( restore_basis_output,arraysize,hami%numberOfVariation )
    read(fd) restore_basis_output(1)%memory(:)
    !!-------
    read(fd) arraysize
    if( allocated(restore_invariant_matrix_spectrum) ) &
         deallocate(restore_invariant_matrix_spectrum)
    allocate( restore_invariant_matrix_spectrum(arraysize) )
    read(fd) restore_invariant_matrix_spectrum(:)
    !!-------
    read(fd) arraysize1,arraysize2,arraysize3
    if( allocated(invariant_matrix_spectrum_initial) ) &
         deallocate(invariant_matrix_spectrum_initial)
    allocate( invariant_matrix_spectrum_initial &
         ( arraysize1, arraysize2, arraysize3 ) )
    read(fd) invariant_matrix_spectrum_initial(:,:,:)
    !!-------
    read(fd) arraysize
    if( allocated(restoreEigenvalue) ) &
         deallocate(restoreEigenvalue)
    allocate( restoreEigenvalue(arraysize) )
    read(fd) restoreEigenvalue(:)
    !!-------
    read(fd) arraysize
    if( allocated(restoreEigenVector) ) &
         deallocate(restoreEigenVector)
    allocate( restoreEigenVector(arraysize) )
    read(fd) restoreEigenVector(:)
    !!-------
    read(fd) arraysize
    if( allocated(reducedDensityMatrix) ) &
         deallocate(reducedDensityMatrix)
    allocate( reducedDensityMatrix(arraysize) )
    read(fd) reducedDensityMatrix(:)
    !!-------
    read(fd) arraysize
    call allocateType( fullSubspaceInfo,arraysize,hami%numberOfVariation )
    read(fd) fullSubspaceInfo(1)%memory(:)
    !!-------
    read(fd) arraysize1, arraysize2
    if( allocated(mapFullSubspaceInfo) ) &
         deallocate(mapFullSubspaceInfo)
    allocate( mapFullSubspaceInfo(arraysize1:arraysize2) )
    read(fd) mapFullSubspaceInfo(:)
    !!-------
    read(fd) arraysize
    read(fd) binnedSpikePositive(:,:)
    !!-------
    read(fd) arraysize
    read(fd) binnedSpikeNegative(:,:)
    !!-------
    close(fd)
  end subroutine loadRestart
end module restart

!********************************************************************
!calculation of thermodynamical properties
!unit of Cn, Sn ::kb
!unit of Xn : kbT/(g_\muB)^2
!all thermodynamical calculation shuld be done before truncation
! %TODO if FLAG_NOIMP set to .true., these properties for free electron chain
!To calculate impurity part (Simp, Ximp etc...), we need to subtract
!the result of FLAG_NOIMP from the result of normal calculation
!******************************************************************
module thermo
  use hamiltonian
  use generall
  use parallel
  use system
  use io_param

  implicit none

  double precision :: betaE

  !partition function
  double precision ::Zn
  !for calculate heat capacity
  double precision :: En
  double precision :: E2n
  !heat capacity
  double precision ::Cn    !<H2>-<H>2
  !entropy
  double precision :: Sn     !E-F

  !for magnetic susceptivility
  double precision ::Sz2n       !<Sz2>
  double precision ::Szn !<Sz>2
  !magnetic susceptivility
  double precision :: Xn

contains
  !calculate the thermodynamical properties in the present iteration
  subroutine calcThermo
    integer ::i,j
    double precision ::expo
    integer ::Sztot
    double precision ::rSztot

    call startCount("calcThermo")

    !initialize
    Zn=0.0d0
    En=0.0d0
    E2n=0.0d0
    Cn=0.0d0
    Sn=0.0d0
    Sz2n=0.0d0
    Szn=0.0d0
    Xn=0.0d0

    do i=1, hami%numberOfBasis
       betaE=betabar*eigenvalue(i)
       expo=exp(-betaE)
       Zn=Zn+expo
       En=En+betaE*expo
       !calculate toatl Sz
       Sztot=0
       do j=hami%conservationBlock(2,1), hami%conservationBlock(2,2)
          Sztot=Sztot+basis_output(i)%basis(j)
       end do
       !Sz=2*Sz(actual)
       rSztot=dble(Sztot)*0.5d0

       Sz2n=Sz2n+expo*rSztot*rSztot
       Szn=Szn+expo*rSztot
    end do

    !use different type of definition of heat capacity based on fluctuation-dissapation theorem
    do i=1, hami%numberOfBasis
       betaE=betabar*eigenvalue(i)
       expo=exp(-betaE)
       Cn=Cn+(betaE-En/Zn)*(betaE-En/Zn)*expo
    end do

    Cn=Cn/Zn
    Sn=En/Zn+log(Zn)
    Xn=Sz2n/Zn-(Szn/Zn)*(Szn/Zn)

    if (my_rank .eq. 0) then
       print*, "Cn=", Cn
    end if

    call stopCount
  end subroutine calcThermo

end module thermo
