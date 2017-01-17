subroutine diagonalization(iteration,ismysubspace)
!!$subroutine diagonalization(iteration,loadParRank)
  use omp_lib
  use generall
  use parallel
  use hamiltonian, only: hami ! in
  use hamiltonian, only: operation ! in
  use hamiltonian, only: subspaceInfo ! in
  use hamiltonian, only: basis_input ! in
  use hamiltonian, only: basis_output ! out
  use hamiltonian, only: eigenvalue ! out
  use hamiltonian, only: eigenvector ! out
  use hamiltonian, only: numberOfSubspace ! in
  use hamiltonian, only: groundStateEnergy ! out
  use hamiltonian, only: sizeOfEigenvector ! out
  use hamiltonian, only: chain_diagonal ! in
  use hamiltonian, only: chain_nondiagonal ! in
  use hamiltonian, only: chain_BCS ! in
  use hamiltonian, only: coefficient_diagonalization_type ! in
  use hamiltonian, only: coefficient_diagonalization ! in
  use hamiltonian, only: coefficient_BCS !in
  use hamiltonian, only: invariant_matrix ! in
  use hamiltonian, only: allocateType
  use system, only: param ! in

  implicit none
  include 'mpif.h'

  integer, intent(in) :: iteration
  logical, intent(in) :: ismysubspace(numberOfSubspace,0:numberOfProcess-1)
!!$  integer, intent(in) :: loadParRank(0:numberOfProcess-1,1:2)

  integer :: j, p, isub, isubl, isubr
  integer :: tmpsizeOfEigenvector, dimSubspace
  double precision, allocatable :: subspaceMatrix(:,:)
  integer :: start_input, start_output, count_eigenvector
  integer :: iinput, jinput, ichain, ioperation, joperation, ireference, jreference, icoef
  integer :: invariant_matrix_type, channel, channelBCS
  double precision :: coefficientNondiagonal
  logical :: flag_coef, flag_BCS_coef
  double precision,allocatable::subspaceEigen(:),work(:)
  external :: dsyev
  integer :: ierr


  integer, allocatable :: tmp(:,:)
  character(20)::Numitr
  integer :: ibl, ibr

  ! MPI variables for eigenvalues
  double precision, allocatable :: mpi_eigenvalue_buffer(:)
  integer, allocatable          :: mpi_eigenvalue_displs(:)
  integer, allocatable          :: omp_eigenvalue_displs(:,:)
  integer, allocatable          :: mpi_eigenvalue_counts(:)
  integer                       :: mpi_eigenvalue_displ
  integer                       :: omp_eigenvalue_displ

  ! MPI variables for eigenvectors
  double precision, allocatable :: mpi_eigenvector_buffer(:)
  integer, allocatable          :: mpi_eigenvector_displs(:)
  integer, allocatable          :: omp_eigenvector_displs(:,:)
  integer, allocatable          :: mpi_eigenvector_counts(:)
  integer                       :: mpi_eigenvector_displ
  integer                       :: omp_eigenvector_displ

  call startCount("diagonal")

  !calculate the size of eigenvector and eigen
  tmpsizeOfEigenvector=0
  do isub=1, numberOfSubspace
     tmpsizeOfEigenvector = tmpsizeOfEigenvector &
          + subspaceInfo(isub)%dimension**2
  end do

  !replace the new information to old information
  sizeOfEigenvector=tmpsizeOfEigenvector

  if( allocated(eigenvector) ) deallocate(eigenvector)
  allocate( eigenvector(sizeOfEigenvector) )

  allocate( mpi_eigenvalue_buffer(hami%numberOfBasis) )
  allocate( mpi_eigenvalue_displs(0:numberOfProcess-1) )
  allocate( omp_eigenvalue_displs(numberOfSubspace,0:numberOfProcess-1) )
  allocate( mpi_eigenvalue_counts(0:numberOfProcess-1) )
  mpi_eigenvalue_buffer = 0.0d0

  allocate( mpi_eigenvector_buffer(sizeOfEigenvector) )
  allocate( mpi_eigenvector_displs(0:numberOfProcess-1) )
  allocate( omp_eigenvector_displs(numberOfSubspace,0:numberOfProcess-1) )
  allocate( mpi_eigenvector_counts(0:numberOfProcess-1) )

  do p=0, numberOfProcess-1
     mpi_eigenvalue_counts(p) = 0
     mpi_eigenvector_counts(p) = 0
     do isub=1, numberOfSubspace
        if( .not. ismysubspace(isub,p) ) cycle
        mpi_eigenvalue_counts(p) = mpi_eigenvalue_counts(p) &
             + subspaceInfo(isub)%dimension
        mpi_eigenvector_counts(p) = mpi_eigenvector_counts(p) &
             + subspaceInfo(isub)%dimension**2
     end do

     if( p==0 ) then
        mpi_eigenvalue_displs(p) = 0
        mpi_eigenvector_displs(p) = 0
     else
        mpi_eigenvalue_displs(p) = mpi_eigenvalue_displs(p-1) &
             + mpi_eigenvalue_counts(p-1)
        mpi_eigenvector_displs(p) = mpi_eigenvector_displs(p-1) &
             + mpi_eigenvector_counts(p-1)
     end if
  end do

  do p=0, numberOfProcess-1
     omp_eigenvalue_displ  = mpi_eigenvalue_displs(p)
     omp_eigenvector_displ  = mpi_eigenvector_displs(p)
     do isub=1, numberOfSubspace
       if( .not. ismysubspace(isub,p) ) cycle
  
       omp_eigenvalue_displs(isub,p) = omp_eigenvalue_displ
       omp_eigenvector_displs(isub,p) = omp_eigenvector_displ

       omp_eigenvalue_displ = omp_eigenvalue_displ &
         + subspaceInfo(isub)%dimension
       omp_eigenvector_displ = omp_eigenvector_displ &
         + subspaceInfo(isub)%dimension**2
    end do
  end do
  


  if( .true. ) then
!!$     mpi_eigenvalue_displ  = mpi_eigenvalue_displs(my_rank)
!!$     mpi_eigenvector_displ = mpi_eigenvector_displs(my_rank)

     !$OMP PARALLEL DO SCHEDULE(dynamic) &
     !$OMP PRIVATE(isub,omp_eigenvalue_displ,omp_eigenvector_displ) &
     !$OMP PRIVATE(dimSubspace,subspaceMatrix,subspaceEigen) &
     !$OMP PRIVATE(start_input,start_output,count_eigenvector) &
     !$OMP PRIVATE(iinput,ioperation,ireference,ichain) &
     !$OMP PRIVATE(jinput,joperation,jreference,flag_coef,icoef) &
     !$OMP PRIVATE(invariant_matrix_type,channel,coefficientNondiagonal) &
     !$OMP PRIVATE(work,isubl,isubr,j,ierr)
     do isub=1, numberOfSubspace
        if( .not. ismysubspace(isub,my_rank) ) cycle

        omp_eigenvalue_displ = omp_eigenvalue_displs(isub,my_rank)
        omp_eigenvector_displ = omp_eigenvector_displs(isub,my_rank)

!!$  if(loadParRank(my_rank,1) .ne. 0 .and. loadParRank(my_rank,2) .ne. 0) then
!!$     do i=loadParRank(my_rank,1), loadParRank(my_rank,2)
        dimSubspace=subspaceInfo(isub)%dimension
        allocate( subspaceMatrix(dimSubspace,dimSubspace) )
        allocate( subspaceEigen(dimSubspace) )

        subspaceMatrix=0.0d0
        start_input=subspaceInfo(isub)%start_input
        start_output=start_input
        count_eigenvector=subspaceInfo(isub)%count_eigenvector

        do iinput=1, dimSubspace
           ioperation=basis_input(start_input+iinput-1)%operation
           ireference=basis_input(start_input+iinput-1)%reference
           !diagonal
           subspaceMatrix(iinput,iinput) = eigenvalue(ireference)*dsqrt(param%Lambda)
           do ichain=1, hami%numberOfChain
              subspaceMatrix(iinput,iinput) = subspaceMatrix(iinput,iinput) &
                   + operation(ioperation, hami%numberOfVariation+ichain) &
                   * chain_diagonal(iteration+1, ichain)
           end do

           if (dimSubspace > 1) then
              !nondiagonal
              !upper triangle
              inner: do jinput=iinput+1, dimSubspace
                 !Here I have to redifine ioperation and ireference, since several case, I switch in the following section
                 ioperation=basis_input(start_input+iinput-1)%operation
                 ireference=basis_input(start_input+iinput-1)%reference
                 joperation=basis_input(start_input+jinput-1)%operation
                 jreference=basis_input(start_input+jinput-1)%reference

                 !making ipoeration < joperation, because only such processes are written in the coefficient file
                 ! due to the Hermite Hamiltonian (tnfn^dagger fn+1 +tnfn+1^dagger fn), inverse process always exist
!                 if (ioperation > joperation) then
!                    ioperation=basis_input(start_input+jinput-1)%operation
!                    ireference=basis_input(start_input+jinput-1)%reference
!                    joperation=basis_input(start_input+iinput-1)%operation
!                    jreference=basis_input(start_input+iinput-1)%reference
!                 end if

                 flag_coef=.false.
                 do icoef=1, hami%numberOfCoefficient
                    if(ioperation .eq. coefficient_diagonalization_type(icoef,1) &
                         .and. joperation .eq. coefficient_diagonalization_type(icoef,2)) then
                       flag_coef=.true.
                       invariant_matrix_type=coefficient_diagonalization_type(icoef,3)
                       channel=coefficient_diagonalization_type(icoef,4)
                       coefficientNondiagonal=coefficient_diagonalization(icoef)
                    else if (ioperation .eq. coefficient_diagonalization_type(icoef,2) &
                         .and. joperation .eq. coefficient_diagonalization_type(icoef,1)) then
                        flag_coef=.true.
                        invariant_matrix_type=coefficient_diagonalization_type(icoef,3)
                        channel=coefficient_diagonalization_type(icoef,4)
                        coefficientNondiagonal=coefficient_diagonalization(icoef)

                        !switch the references
                        ireference=basis_input(start_input+jinput-1)%reference
                        jreference=basis_input(start_input+iinput-1)%reference
                    end if
                 end do

                 !BCS gap
                 if (hami%flag_BCS) then
                 do icoef=1, hami%numberOfBCSCoefficient
                    if ((ioperation .eq. coefficient_BCS(icoef,1)) .and. (joperation .eq. coefficient_BCS(icoef,2)) .and. (ireference .eq. jreference)) then
                        frag_BCS_coef=.true.
                        channelBCS=coefficient_BCS(icoef,3)

                    end if
                 end do
                 end if

                 if((.not.flag_coef) .and. (.not. flag_BCS_coef)) cycle inner

                 !normal lead
                 if (flag_coef) then
                 subspaceMatrix(iinput,jinput) &
                      = chain_nondiagonal(iteration, channel) &
                      * coefficientNondiagonal &
                      * invariant_matrix(jreference,ireference,invariant_matrix_type)

                 subspaceMatrix(jinput,iinput) &
                      = subspaceMatrix(iinput,jinput)
                 end if

                 !BCS paring
                 if (hami%flag_BCS) then
                    if (flag_BCS_coef) then
                    subspaceMatrix(iinput, jinput) =chain_bcs(iteration+1, channelBCS)
                    subspaceMatrix(jinput, iinput) =subspaceMatrix(iinput, jinput)
                    end if
                 end if

!!$                 !debug for general QSz
!!$                 if (iteration==0 .and. isub==11) then
!!$                    write(buglog,*) "invariantmatrix element", &
!!$                         invariant_matrix(jreference,ireference,invariant_matrix_type)
!!$                    write(buglog,*) "chain_nondiagonal", &
!!$                         chain_nondiagonal(iteration, channel)
!!$                    write(buglog,*) "coefficient_diagonalization", &
!!$                         coefficientNondiagonal
!!$                    write(buglog,*) "matrix index i,j=", &
!!$                         iinput, jinput
!!$                    write(buglog,*) "matrix element=", &
!!$                         subspaceMatrix(iinput,jinput)
!!$                 end if
              end do inner
           end if
        end do ! iinput

        !TODO diagonalization
        allocate( work(10*dimSubspace) )
        call dsyev( "V", "U", dimSubspace, subspaceMatrix, dimSubspace, &
             subspaceEigen, work, size(work), ierr )
        deallocate(work)

        !for reduce the numerical error and artificial lift of degeneracy,
        !here I chop the numerical error term in subspaceMatrix
        do isubl=1, dimSubspace
           do isubr=1, dimSubspace
              if (abs(subspaceMatrix(isubl,isubr)) &
                   < eigenvectorThresould) then
                 subspaceMatrix(isubl,isubr) = 0.0d0
              end if
           end do ! isubr
        end do ! isubl

        !filling to container
        do j=1, dimSubspace
          mpi_eigenvalue_buffer(omp_eigenvalue_displ+j) &
                = subspaceEigen(j)

!!$           basis_output(start_output+j-1)%basis(:)=subspaceInfo(i)%basis(:)
!!$           basis_output(start_output+j-1)%reference = j
!!$           basis_output(start_output+j-1)%variation = 0
!!$           basis_output(start_output+j-1)%operation = 0

!!$           eigenvalue(start_output+j-1) = subspaceEigen(j)
        end do ! j
!!$        mpi_eigenvalue_displ = mpi_eigenvalue_displ + dimSubspace

        j=1
        do isubr=1,dimSubspace
           do isubl=1, dimSubspace
              mpi_eigenvector_buffer(omp_eigenvector_displ+j) &
                   = subspaceMatrix(isubl,isubr)
              j=j+1
!!$              eigenvector(count_eigenvector)=subspaceMatrix(isubl,isubr)
!!$              count_eigenvector=count_eigenvector+1
           end do ! isubl
        end do ! isubr
!!$        mpi_eigenvector_displ = mpi_eigenvector_displ + dimSubspace**2

        deallocate(subspaceEigen)
        deallocate(subspaceMatrix)
     end do ! isub
     !$OMP END PARALLEL DO
  end if



  !! allgather eigenvalues
  call MPI_ALLGATHERV( &
       mpi_eigenvalue_buffer(mpi_eigenvalue_displs(my_rank)+1), &
       mpi_eigenvalue_counts(my_rank), &
       MPI_DOUBLE_PRECISION, &
       mpi_eigenvalue_buffer, & 
       mpi_eigenvalue_counts(:), &
       mpi_eigenvalue_displs(:), &
       MPI_DOUBLE_PRECISION, &
       MPI_COMM_WORLD, ierr )

  if( allocated(eigenvalue) ) deallocate(eigenvalue)
  allocate( eigenvalue(hami%numberOfBasis) )
  do p=0, numberOfProcess-1
     mpi_eigenvalue_displ = mpi_eigenvalue_displs(p)
     do isub=1, numberOfSubspace
        if( .not. ismysubspace(isub,p) ) cycle

        dimSubspace  = subspaceInfo(isub)%dimension
        start_output = subspaceInfo(isub)%start_input
        do j=1, dimSubspace
           eigenvalue(start_output+j-1) &
                = mpi_eigenvalue_buffer(mpi_eigenvalue_displ+j)
        end do
        mpi_eigenvalue_displ = mpi_eigenvalue_displ + dimSubspace
     end do
  end do


  deallocate(mpi_eigenvalue_displs)
  deallocate(mpi_eigenvalue_counts)
  deallocate(mpi_eigenvalue_buffer)


  !make the most small eigenvalue to zero
  groundStateEnergy=eigenvalue(1)

  eigenvalue(:) = eigenvalue(:) - minval(eigenvalue(:))
  where( eigenvalue(:) <= numericalThresould )
     eigenvalue(:) = 0.0d0
  end where

  !! allgather eigenvectors
  call MPI_ALLGATHERV( &
       mpi_eigenvector_buffer(mpi_eigenvector_displs(my_rank)+1), &
       mpi_eigenvector_counts(my_rank), &
       MPI_DOUBLE_PRECISION, &
       mpi_eigenvector_buffer, & 
       mpi_eigenvector_counts(:), &
       mpi_eigenvector_displs(:), &
       MPI_DOUBLE_PRECISION, &
       MPI_COMM_WORLD, ierr )

  do p=0, numberOfProcess-1
     mpi_eigenvector_displ = mpi_eigenvector_displs(p)
     do isub=1, numberOfSubspace
        if( .not. ismysubspace(isub,p) ) cycle

        dimSubspace = subspaceInfo(isub)%dimension
        count_eigenvector = subspaceInfo(isub)%count_eigenvector

        j=1
        do isubr=1, dimSubspace
           do isubl=1, dimSubspace
              eigenvector(count_eigenvector+j-1) &
                   = mpi_eigenvector_buffer(mpi_eigenvector_displ+j)
              j=j+1
           end do
        end do
        mpi_eigenvector_displ = mpi_eigenvector_displ + dimSubspace**2
     end do
  end do

  deallocate(mpi_eigenvector_displs)
  deallocate(mpi_eigenvector_counts)
  deallocate(mpi_eigenvector_buffer)


  !! update basis_output
  call allocateType( basis_output, &
       hami%numberOfBasis, hami%numberOfVariation )
  do isub=1, numberOfSubspace
     dimSubspace  = subspaceInfo(isub)%dimension
     start_output = subspaceInfo(isub)%start_input
     do j=1, dimSubspace
        basis_output(start_output+j-1)%basis(:) &
             = subspaceInfo(isub)%basis(:)
        basis_output(start_output+j-1)%reference = j
        basis_output(start_output+j-1)%variation = 0
        basis_output(start_output+j-1)%operation = 0
     end do
  end do

  call stopCount
end subroutine diagonalization
