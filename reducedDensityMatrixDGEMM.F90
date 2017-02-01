!change the algoritm of calculation in reducedDensityMatrix.F90, I will move this file after test
subroutine reducedDensityMatrixCalculation
  use omp_lib
  use generall
  use hamiltonian, only: hami ! in
  use hamiltonian, only: operation ! in
  use system, only: fullSubspaceInfo ! inout
  use system, only: fullBasisInput ! in
  use system, only: fullBasisOutput ! in
  use system, only: elementNumFullSubspace ! in
  use system, only: mapFullSubspaceInfo ! in
  use system, only: reducedDensityMatrix ! inout
  use system, only: pointReducedDensity ! inout
  use system, only: fullEigenvector ! in
  use parallel
  implicit none
  include 'mpif.h'

  integer :: i,j,iS, iSprev, ir, il, ioperation, joperation,iv, ilprev, irprev, ikeepPrev1, ikeepPrev2, idiag, ip
  integer :: ierr
  integer :: presentIteration, previousIteration
  integer :: presentStartSubspace, previousStartSubspace, previousEndSubspace, presentStartOutput
  integer :: presentNumberOfBasis, presentOutput
  integer :: presentIndexLeft, presentIndexRight, presentReferenceLeft, presentReferenceRight
  integer :: previousKeepBasis(hami%numberOfOperation), previousTotalBasis(hami%numberOfOperation)
  integer :: previousPositionReducedDensity(hami%numberOfOperation), previousPositionOutput(hami%numberOfOperation)
  integer :: previousPositionEigenvector(hami%numberOfOperation), previousPositionInput(hami%numberOfOperation)
  integer :: derivedOperation(hami%numberOfOperation), flag_find_space(hami%numberOfOperation)
  integer :: totalBasisNumPrevious
  integer :: reducedDensity_reference
  integer, allocatable :: bcast_start(:), bcast_end(:), numelem(:)
  logical :: flag_find_space
  integer, allocatable :: numbasis(:)
  integer :: nsub !number of subspace in this iteration
  integer :: subcount, total, diff
  integer ::loadParRank(0:numberOfProcess-1,1:2)
  integer,allocatable :: presentVariation(:), previousVariation(:) ! delete this variables, previousVariationL(:)
  double precision, allocatable :: reducedDensitySubspace(:,:)
  double precision, allocatable :: reducedDensityPrevious(:,:), reducedDensityPreviousSubSpace(;,;)
  double precision :: diag, diagtotal
  double precision:: element
  integer :: ilderiv, irderiv
  integer :: IC

  double precision, allocatable :: vectorL(:,:), vectorR(:,:), coefMatrix(:,:,:)
  external :: dgemm
  double precision, parameter :: alpha=1.0, beta=0.0
  integer:: rowA, columnA, rowB, columnB, rowC, columnC
  !integer ::loadmax,loadmin

  ! print*, "numberOfVariation=", numberOfVariation

  call startCount("redDenMat")

  if (my_rank .eq. 0) print*, "start the calculation of reduced density matrix"

  allocate( previousVariation(hami%numberOfVariation) )
  !allocate(previousVariationL(1:hami%numberOfVariation))

  !outer loop for iteration
  do i=1,hami%numberOfIteration-hami%startTrancation
     call startCount("redDenMat:outloop")
     diag=0.0d0

     presentIteration=hami%numberOfIteration-i
     previousIteration=hami%numberOfIteration-i+1
     if(my_rank .eq. 0) then
        print*, "iteration=", presentIteration
     end if

     presentStartSubspace=mapFullSubspaceInfo(presentIteration)
     previousStartSubspace=mapFullSubspaceInfo(previousIteration)
     presentStartOutput=fullSubspaceInfo(presentStartSubspace)%start_output

     !counding the number of basis in each subspace in this iteration
     !this procedure enable the loadbalance and update the starting point of reduced density matrix
     !in each subspace
     !fullSubspaceInfo
     !numberOfVariation+1 : the degeneracy of this subspace before truncation
     !numberOfVariation+2 : degeneracy of this subspace after truncation
     !numberOfVariation+3 : the starting point of this subspace in fullEigenvector
     !numberOfVariation+4 : the starting point of this subspace in fullBasisInput
     !numberOfVariation+5 : the starting point of this subspace in fullBasisOutput (& fullEigenvalue)
     !numberOfVariation+6 : the starting point of this subspace in reducedDensityMatrix, filled in the calculation of reduced density matrix part
     nsub=previousStartSubspace-presentStartSubspace
     if (my_rank .eq. 0) then
        print*, "nsub=", nsub
     end if

     allocate( numbasis(nsub) )
     subcount=1
     !fillup the final element of fullSubspaceInfo
     do j=presentStartSubspace, previousStartSubspace-1
        numbasis(subcount) = fullSubspaceInfo(j)%dimension_after
        fullSubspaceInfo(j)%start_matrix = pointReducedDensity
        pointReducedDensity=pointReducedDensity+fullSubspaceInfo(j)%dimension_after &
             * fullSubspaceInfo(j)%dimension_after
        subcount=subcount+1
     end do

     !set the load barance
     total=0
     do j=1, nsub
        total=total+numbasis(j)
     end do

     call loadbalance1DSimple(loadParRank,nsub)

     call stopCount
     call startCount("redDenMat:inloop")


     !parallel part
     !$OMP PARALLEL DO &
     !$OMP SCHEDULE(dynamic) &
     !$OMP PRIVATE(iS, iSprev, ir, il, ioperation, joperation, iv, ilprev, irprev)&
     !$OMP PRIVATE(ikeepPrev1, ikeepPrev2, idiag, previousEndSubspace, derivedOperation)&
     !$OMP PRIVATE(presentNumberOfBasis, presentOutput, presentIndexLeft, presentIndexRight)&
     !$OMP PRIVATE(presentReferenceLeft, presentReferenceRight, totalBasisNumPrevious )&
     !$OMP PRIVATE(previousKeepBasis, previousTotalBasis, previousPositionReducedDensity)&
     !$OMP PRIVATE(previousPositionOutput, previousPositionInput, previousPositionEigenvector)&
     !$OMP PRIVATE(reducedDensity_reference, flag_find_space,diff, presentVariation)&
     !$OMP PRIVATE(previousVariation, reducedDensitySubspace, element, ilderiv, irderiv)&
     !$OMP PRIVATE(reducedDensityPrevious, reducedDensityPreviousSubSpace)&
     !$OMP REDUCTION(+: diag)
     do iS=presentStartSubspace+loadParRank(my_rank,1)-1, presentStartSubspace+loadParRank(my_rank,2)-1
        !$ IC=omp_get_thread_num()
        !print*, "load::",iS, "thread::", IC
        !print*, "iS=", iS
        !TODO reduced density matrix calculation
        allocate( presentVariation(hami%numberOfVariation) )
        do iv=1, hami%numberOfVariation
           presentVariation(iv)=fullSubspaceInfo(iS)%basis(iv)
        end do

        presentNumberOfBasis=fullSubspaceInfo(iS)%dimension_after
        presentOutput=fullSubspaceInfo(iS)%start_output

        allocate( reducedDensitySubspace(presentNumberOfBasis,presentNumberOfBasis) )
        reducedDensitySubspace=0.0d0


             !initialize
             flag_find_space=.false.
             previousTotalBasis=0
             previousKeepBasis=0
             previousPositionEigenvector=0
             previousPositionInput=0
             previousPositionOutput=0
             previousPositionReducedDensity=0
             derivedOperation=0


             !find related subspace in the pervious iteration and construct vectorL
             !loop_operation: do ioperation=1, numberOfOperation
              do ioperation =1, hami%numberOfOperation

                 do iv=1, hami%numberOfVariation
                    !caution :: present iteration=N previous iteration=N+1,
                    !the procedure of calculation of reduced density matrix is backword of the normal NRG calculation
                    previousVariation(iv)=presentVariation(iv)-operation(ioperation , iv)
                 end do

                 flag_find_space=.false.

                 !find corresponding subspace in previousIteration
                 if(previousIteration .eq. hami%numberOfIteration) then
                    !for last iteration (first iteration for reduced density matrix calc), the end point is the end of container
                    previousEndSubspace=elementNumFullSubspace
                 else
                    previousEndSubspace=mapFullSubspaceInfo(previousIteration+1)-1
                 end if

                 !serch corresponding subspace
                 do iSprev=mapFullSubspaceInfo(previousIteration), previousEndSubspace
                    diff=0
                    do iv=1, hami%numberOfVariation
                       diff=diff+abs(previousVariation(iv)-fullSubspaceInfo(iSprev)%basis(iv))
                    end do

                    if (diff .eq. 0) then
                       flag_find_space(ioperation)=.true.
                       previousTotalBasis(ioperation)=fullSubspaceInfo(iSprev)%dimension_before
                       previousKeepBasis(ioperation)=fullSubspaceInfo(iSprev)%dimension_after
                       previousPositionEigenvector(ioperation)=fullSubspaceInfo(iSprev)%start_eigen
                       previousPositionInput(ioperation)=fullSubspaceInfo(iSprev)%start_input
                       previousPositionOutput(ioperation)=fullSubspaceInfo(iSprev)%start_output
                       previousPositionReducedDensity(ioperation)=fullSubspaceInfo(iSprev)%start_matrix
                       derivedOperation(ioperation)=ioperation
                       exit
                    end if
                 end do !do iSprev, to find the corresponding subspace
            end do ! operation

            !count the number of total basis derived from this subspace in n+1 shell
            totalBasisNumPrevious=sum(previousKeepBasis)

            do ioperation=1, hami%numberOfOperation
            if (flag_find_space(ioperation)) then
                allocate(reducedDensityPreviousSubSpace(previous))
                do il=1, presentNumberOfBasis
                presentIndex=fullBasisOutput(presentOutput+il-1)%reference
                    do ilprev=previousPositionInput,previousPositionInput+previousTotalBasis-1
                       if((fullBasisInput(ilprev)%variation .eq. presentIndexLeft) &
                            .and. (fullBasisInput(ilprev)%operation .eq. ioperation) ) then
                          !the position of derived state from il-th basis
                          ilderiv=ilprev
                          exit
                       end if
                    end do




                end do

            end if
            end do


        !double loop for each basis

        !for calculate matrix element reducedDensitySubspace(il,ir)
        do il=1, presentNumberOfBasis
           !for basis_output
           !element placed at numberOfVariation+1 : variation : index of basis with same conserved quantity
           presentIndexLeft=fullBasisOutput(presentOutput+il-1)%reference

           !add routine to calculate the reference correspond to reference in fullbasisInput
           presentReferenceLeft=presentOutput-presentStartOutput+il
           do ir=1,presentNumberOfBasis

              presentIndexRight=fullBasisOutput(presentOutput+ir-1)%reference
              !add routine to calculate the reference correspond to reference in fullbasisInput
              presentReferenceRight=presentOutput-presentStartOutput+ir

              !initialize the matrix element : rho_(il, ir)
              element=0.0d0

              !loop_operation: do ioperation=1, numberOfOperation
              do ioperation =1, hami%numberOfOperation

                 do iv=1, hami%numberOfVariation
                    !caution :: present iteration=N previous iteration=N+1,
                    !the procedure of calculation of reduced density matrix is backword of the normal NRG calculation
                    previousVariation(iv)=presentVariation(iv)-operation(ioperation , iv)
                 end do

                 flag_find_space=.false.

                 !find corresponding subspace in previousIteration
                 if(previousIteration .eq. hami%numberOfIteration) then
                    !for last iteration (first iteration for reduced density matrix calc), the end point is the end of container
                    previousEndSubspace=elementNumFullSubspace
                 else
                    previousEndSubspace=mapFullSubspaceInfo(previousIteration+1)-1
                 end if

                 !serch corresponding subspace
                 do iSprev=mapFullSubspaceInfo(previousIteration), previousEndSubspace
                    diff=0
                    do iv=1, hami%numberOfVariation
                       diff=diff+abs(previousVariation(iv)-fullSubspaceInfo(iSprev)%basis(iv))
                    end do

                    if (diff .eq. 0) then
                       flag_find_space=.true.
                       previousTotalBasis=fullSubspaceInfo(iSprev)%dimension_before
                       previousKeepBasis=fullSubspaceInfo(iSprev)%dimension_after
                       previousPositionEigenvector=fullSubspaceInfo(iSprev)%start_eigen
                       previousPositionInput=fullSubspaceInfo(iSprev)%start_input
                       previousPositionOutput=fullSubspaceInfo(iSprev)%start_output
                       previousPositionReducedDensity=fullSubspaceInfo(iSprev)%start_matrix
                       exit
                    end if
                 end do !do iSprev, to find the corresponding subspace


                 if (flag_find_space) then

                    !for basis_input
                    !element placed at numberOfVariation+1 : reference : index to define from which eigenstate of last iteration
                    !element placed at numberOfVariation+2 : variation : index of basis with same conserved quantity
                    !element placed at numberOfVariation+3 : operation :what kind of operation is done for creating this basis from the eigenstate of last iteration (add f_{ni\sigma}^dagger etc)

                    !determine the state derived from il-th basis with this operation
                    do ilprev=previousPositionInput,previousPositionInput+previousTotalBasis-1
                       if((fullBasisInput(ilprev)%variation .eq. presentIndexLeft) &
                            .and. (fullBasisInput(ilprev)%operation .eq. ioperation) ) then
                          !the position of derived state from il-th basis
                          ilderiv=ilprev
                          exit
                       end if
                    end do

                    !loop for the making derivative state with ir-th basis
                    joperation=ioperation
                    do irprev=previousPositionInput,previousPositionInput+previousTotalBasis-1
                       if((fullBasisInput(irprev)%variation .eq. presentIndexRight) &
                            .and. (fullBasisInput(irprev)%operation .eq. joperation) ) then
                          irderiv=irprev
                          exit
                       end if
                    end do

                    !summing up the contribution from each basis_input by using the product of eigenvector and reduced density

                    do ikeepPrev1=1, previousKeepBasis
                       do ikeepPrev2=1, previousKeepBasis
                          element=element+fullEigenvector(previousPositionEigenvector &
                               + (ikeepPrev1-1)*previousTotalBasis+ilderiv-previousPositionInput) &
                               * fullEigenvector(previousPositionEigenvector &
                               + (ikeepPrev2-1)*previousTotalBasis+irderiv-previousPositionInput) &
                               * reducedDensityMatrix(previousPositionReducedDensity &
                               + (ikeepPrev1-1)*previousKeepBasis+ikeepPrev2-1)
                       end do
                    end do
                 end if ! if for find subspace derived from il-th basis
              end do !making derivative by operation ioperation

              reducedDensitySubspace(il,ir) = element
           end do
        end do

        !before exiting  subspace
        !for check the normalization
        do idiag=1,presentNumberOfBasis
           diag=diag+reducedDensitySubspace(idiag,idiag)
        end do

        !fill reducedDensitMatrix
        reducedDensity_reference=fullSubspaceInfo(iS)%start_matrix
        do il=1,presentNumberOfBasis
           do ir=1,presentNumberOfBasis
              reducedDensityMatrix(reducedDensity_reference)=reducedDensitySubspace(il,ir)
              reducedDensity_reference=reducedDensity_reference+1
           end do
        end do

        deallocate(reducedDensitySubspace)
        deallocate(presentVariation)
     end do
     !$OMP END PARALLEL DO

     call stopCount

     if (my_rank .eq. 0) then
        print*, "summation of diagonal element in each rank=", diag
     end if

     call startCount("redDenMat:Allgather")

     allocate( bcast_start(0:numberOfProcess-1) )
     allocate( bcast_end(0:numberOfProcess-1) )
     allocate( numelem(0:numberOfProcess-1) )

     !TODO broadcast
     !set the start and end point of broadcast
     do ip=0, numberOfProcess-1
        bcast_start(ip)=fullSubspaceInfo(presentStartSubspace+loadParRank(ip,1)-1)%start_matrix
        bcast_end(ip)=fullSubspaceInfo(presentStartSubspace+loadParRank(ip,2)-1)%start_matrix &
             +fullSubspaceInfo(presentStartSubspace+loadParRank(ip,2)-1)%dimension_after**2-1

        !20130409 fix the count of numelem
        numelem(ip)=bcast_end(ip)-bcast_start(ip)+1
!!$        !print*, "bcast_start, bcast_end,numelem", bcast_start, bcast_end,numelem
!!$        call MPI_BCAST(reducedDensityMatrix(bcast_start(ip)), numelem(ip), MPI_DOUBLE_PRECISION, ip, MPI_COMM_WORLD, ierr)
     end do

     call MPI_ALLGATHERV( &
          reducedDensityMatrix(bcast_start(my_rank)), &
          numelem(my_rank), &
          MPI_DOUBLE_PRECISION, &
          reducedDensityMatrix, &
          numelem(:), &
          bcast_start(:)-1, &
          MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr )

     deallocate( bcast_start )
     deallocate( bcast_end )
     deallocate( numelem )
     call stopCount

     call startCount("redDenMat:Reduce")
     call MPI_REDUCE( diag, diagtotal, &
          1, MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD, ierr )

     if(my_rank .eq. 0) print*, "total sum of diagonal element=", diagtotal

     deallocate(numbasis)

     call stopCount
  end do

  if (my_rank .eq. 0) print*, "finish :: reduced density matrix"

  call stopCount
end subroutine reducedDensityMatrixCalculation
