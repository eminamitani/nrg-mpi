subroutine invariantMatrixForSpectrumDGEMM
    use omp_lib
    use generall
    use parallel
    use hamiltonian, only: SubspaceInfo_type ! type
    use hamiltonian, only: hami ! in
    use hamiltonian, only: invariant_matrix_spectrum ! inout
    use hamiltonian, only: subspaceInfo ! in
    use hamiltonian, only: basis_output ! in
    use hamiltonian, only: basis_input ! in
    use hamiltonian, only: numberOfSubspace ! in
    use hamiltonian, only: coefficient_invariant_matrix_spectrum ! in
    use hamiltonian, only: eigenvector ! in
    use hamiltonian, only: conservation_difference_spectrum ! in

    implicit none
  include 'mpif.h'

    integer :: loadmin(0:numberOfProcess-1), myloadmin
    integer :: loadmax(0:numberOfProcess-1), myloadmax
    integer :: allload(0:numberOfProcess-1), myload

    integer :: ileft, iright
    type(SubspaceInfo_type) :: subspace_left, subspace_right
    integer :: chargeStart, chargeEnd, spinStart, spinEnd
    integer, allocatable :: difference(:)
    logical, allocatable :: flags(:)
    integer :: imat
    logical :: flag_total
    integer :: countTrue
    integer, allocatable :: countMatrix(:)
    integer :: countIndex
    integer :: isub, ipair
    integer, allocatable :: tmp2(:)
    integer :: rmax_left, rmax_right
    integer :: eigenvec_min_left,eigenvec_min_right
    integer :: basis_min_left, basis_min_right
    integer :: basis_max_left, basis_max_right
    integer :: isubl, isubr
    integer :: variation_left, variation_right
    integer :: operation_left, operation_right
    integer :: reference_left, reference_right
    integer :: diff_var_left, diff_var_right
    integer :: imatrix, startl, startr
    double precision :: coef

    integer, allocatable :: subspacePair(:,:)

    integer, allocatable :: keeped_basis_number(:)
    integer :: degeneracy_truncation, startOfBasisOutput

    logical :: flag_parallel !flag for parallel or not
    integer :: ierr
    integer :: p, m, i,j

    double precision, allocatable :: invariant_matrix_spectrum_new(:,:,:) !the renew value

    double precision, allocatable :: vectorL(:,:), vectorR(:,:), coefMatrix(:,:,:)
    !intermediate matrix for dgemm
    double precision, allocatable :: tmp1(:,:)
    !invariant matrix for subspace
    double precision, allocatable :: invariant_matrix_sub(:,:)
    external :: dgemm
    double precision, parameter :: alpha=1.0, beta=0.0
    integer:: rowA, columnA, rowB, columnB, rowC, columnC

    call startCount("invMatSpec")
    call startCount("invMatSpec:omp")
    !if the numberOfBasis < numberOfProcess, avoid to paralell
    allocate(keeped_basis_number(numberOfSubspace))

    do isub=1,  numberOfSubspace
        degeneracy_truncation=0
        startOfBasisOutput=1

        do j=1, hami%numberOfBasis
            if( all(subspaceInfo(isub)%basis(:)==basis_output(j)%basis(:)) ) then
                if(startOfBasisOutput > j) then
                    startOfBasisOutput=j
                end if
                degeneracy_truncation=degeneracy_truncation+1
            end if
        end do
        keeped_basis_number(isub)=degeneracy_truncation

    end do

    !print*,keeped_basis_number

    !prepare the subspace pair set
    if( allocated(subspacePair) ) deallocate(subspacePair)
    allocate (subspacePair(numberOfSubspace*numberOfSubspace,2))

    do i=1, numberOfSubspace
        do j=1,numberOfSubspace
            subspacePair((i-1)*numberOfSubspace+j,1)=i
            subspacePair((i-1)*numberOfSubspace+j,2)=j
        end do
    end do

       !loadbalancing

    if (numberOfSubspace*numberOfSubspace < numberOfProcess*4 ) then
        if (my_rank .eq. 0) print*, "invariantMatrix:: too small system size, avoid parallel process"
        !do calculation for all matrix in each rank
        myloadmin = 1
        myloadmax = numberOfSubspace*numberOfSubspace
        myload    = myloadmax - myloadmin + 1
        flag_parallel= .false.
    else
        call loadbalance3( loadmin, loadmax, allload, numberOfSubspace*numberOfSubspace )
        myloadmin = loadmin(my_rank)
        myloadmax = loadmax(my_rank)
        myload    = allload(my_rank)
        flag_parallel= .true.
    end if

    allocate( invariant_matrix_spectrum_new &
        ( hami%numberOfBasis, hami%numberOfBasis, hami%numberOfMatrix) )

    invariant_matrix_spectrum_new = 0.0d0


    !$OMP PARALLEL DO &
    !$OMP SCHEDULE(dynamic) &
    !$OMP PRIVATE(ileft,iright, ipair) &
    !$OMP PRIVATE(subspace_left,subspace_right) &
    !$OMP PRIVATE(chargeStart,chargeEnd,spinStart,spinEnd) &
    !$OMP PRIVATE(difference,flags,imat) &
    !$OMP PRIVATE(flag_total,countTrue) &
    !$OMP PRIVATE(countMatrix,countIndex) &
    !$OMP PRIVATE(isub) &
    !$OMP PRIVATE(rmax_left,rmax_right) &
    !$OMP PRIVATE(eigenvec_min_left,eigenvec_min_right) &
    !$OMP PRIVATE(basis_min_left,basis_min_right) &
    !$OMP PRIVATE(basis_max_left,basis_max_right) &
    !$OMP PRIVATE(isubl,isubr) &
    !$OMP PRIVATE(variation_left,variation_right) &
    !$OMP PRIVATE(operation_left,operation_right) &
    !$OMP PRIVATE(reference_left,reference_right) &
    !$OMP PRIVATE(diff_var_left,diff_var_right) &
    !$OMP PRIVATE(imatrix,coef)&
    !$OMP PRIVATE(vectorL, vectorR, coefMatrix, tmp1, invariant_matrix_sub)&
    !$OMP PRIVATE(startl,startr)&
    !$OMP PRIVATE(rowA, columnA, rowB, columnB, rowC, columnC)

    loop_pair:do ipair=myloadmin, myloadmax
        ileft=subspacePair(ipair,1)
        iright=subspacePair(ipair,2)
        if ((keeped_basis_number(ileft) .eq. 0) .or. (keeped_basis_number(iright) .eq. 0) ) cycle loop_pair

        subspace_left = subspaceInfo(ileft)
        subspace_right = subspaceInfo(iright)

        
        !calculate <ileft | d^dagger |iright> type  invariant matrix element for single particle excitation spectrum

        !except the physically unmeaning pattern,
        !this depends on the type of model Hamiltonian


        allocate( difference(hami%numberOfVariation) )
        allocate( flags(hami%numberOfMatrix) )

        difference(:) = subspace_left%basis(:)-subspace_right%basis(:)

        do imat=1, hami%numberOfMatrix
            flags(imat) = all( difference(:) == conservation_difference_spectrum(imat,:) )
        end do
        deallocate(difference)

        flag_total=.true.
        countTrue=0
        do imat=1, hami%numberOfMatrix
            flag_total=flag_total.and.flags(imat)
            if (flags(imat)) countTrue=countTrue+1
        end do

        if (countTrue >0) then
            allocate( countMatrix(countTrue) )
            countIndex=1
            do imat=1, hami%numberOfMatrix
                if (flags(imat)) then
                    countMatrix(countIndex)=imat
                    countIndex=countIndex+1
                end if
            end do
        end if

        deallocate(flags)
        if (countTrue == 0) then
            cycle loop_pair
        end if


        rmax_right = subspaceInfo(iright)%dimension
        eigenvec_min_right = subspaceInfo(iright)%count_eigenvector
        basis_min_right = subspaceInfo(iright)%start_input
        basis_max_right = basis_min_right+rmax_right-1


        rmax_left = subspaceInfo(ileft)%dimension
        eigenvec_min_left = subspaceInfo(ileft)%count_eigenvector
        basis_min_left = subspaceInfo(ileft)%start_input
        basis_max_left = basis_min_left+rmax_left-1


        !constructing matrix of eigenvector and coefficient matrix
        !Invariant matrix = (matrix of eigenvector L)^T * coefficient matrix * (matrix of eigenvector R)
        ! Here, vectorL is stored as (matrix of eigenvector L)^T
        ! and vectorR is stored as matrix of eigenvector R


        !vector L transposed
        allocate (vectorL(keeped_basis_number(ileft), rmax_left))
        vectorL=0.0


        allocate (vectorR(rmax_right,keeped_basis_number(iright)))
        vectorR=0.0

        do isubl=1, keeped_basis_number(ileft)
            vectorL(isubl, 1:rmax_left)=eigenvector(eigenvec_min_left+(isubl-1)*rmax_left:eigenvec_min_left+isubl+rmax_left-1)
        end do

        do isubr=1, keeped_basis_number(iright)
            vectorR(1:rmax_right, isubr)=eigenvector(eigenvec_min_right+(isubr-1)*rmax_right:eigenvec_min_right+isubr+rmax_right-1)
        end do


        allocate (coefMatrix(rmax_left,rmax_right,hami%numberOfConductionMatrix))

        do isubl=basis_min_left, basis_max_left

            operation_left=basis_input(isubl)%operation

            loop_sub_r: do isubr=basis_min_right, basis_max_right

                operation_right=basis_input(isubr)%operation

                if (operation_right .ne. operation_left) then
                    cycle loop_sub_r
                end if

                do imat=1, size(countMatrix)
                    imatrix = countMatrix(imat)
                    coefMatrix(isubl-basis_min_left+1,isubr-basis_min_right+1, imatrix) = coefficient_invariant_matrix_spectrum &
                        (imatrix,operation_left)
                end do

            end do loop_sub_r
        end do

        allocate (tmp1(keeped_basis_number(ileft),rmax_right ))


        allocate (invariant_matrix_sub(keeped_basis_number(ileft), keeped_basis_number(iright)))


        do imat=1, size(countMatrix)
            imatrix = countMatrix(imat)
            tmp1=0.0
            invariant_matrix_sub=0.0

            !print*,"matrix kind:", ip
            !print*, coefMatrix(:,:,ip)

            !vectorL * coefMatrix
            rowA=keeped_basis_number(ileft)
            columnA=rmax_left
            rowB=rmax_left
            columnB=rmax_right
            rowC=keeped_basis_number(ileft)
            columnC=rmax_right

            call dgemm('N','N',rowA,columnB,columnA,alpha,vectorL,rowA,&
                coefMatrix(:,:,imatrix),rowB,beta,tmp1,rowC)
            !print*,"tmp1"
            !print*,tmp1

            !tmp1--> A
            rowA=keeped_basis_number(ileft)
            columnA=rmax_right
            !vectorR-->B
            rowB=rmax_right
            columnB=keeped_basis_number(iright)
            !invariant -->C
            rowC=keeped_basis_number(ileft)
            columnC=keeped_basis_number(iright)

            !tmp1*vectorR^T
            call dgemm('N','N',rowA, columnB,columnA,alpha,tmp1,rowA,&
                vectorR,rowB,beta,invariant_matrix_sub,rowC)
            !print*,"invariant"
            !print*, invariant_matrix_sub

            !copy to invariant Matrix
            startl=sum(keeped_basis_number(1:ileft-1))+1
            startr=sum(keeped_basis_number(1:iright-1))+1

            !print*,"startl, startr:", startl, " ", startr

            do isubl=1, keeped_basis_number(ileft)
                do isubr=1, keeped_basis_number(iright)
                    !for memory access, invariant_matrix index is reverted
                    invariant_matrix_spectrum_new(startr+isubr-1, startl+isubl-1, imatrix)=invariant_matrix_sub(isubl,isubr)
                end do
            end do


        end do

        deallocate(vectorL, vectorR, coefMatrix, tmp1, invariant_matrix_sub)

        deallocate(countMatrix)
    end do loop_pair

    !$OMP END PARALLEL DO

    call stopCount
    call startCount("invMatSpec:Allgather")

    if (flag_parallel) then
!        do m=1, hami%numberOfMatrix
!            call MPI_ALLGATHERV( &
!                invariant_matrix_spectrum_new(1,myloadmin,m), &
!                allload(my_rank)*hami%numberOfBasis, &
!                MPI_DOUBLE_PRECISION, &
!                invariant_matrix_spectrum_new(1,1,m), &
!                allload(:)*hami%numberOfBasis, &
!                (loadmin(:)-1)*hami%numberOfBasis, &
!                MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr )
!        end do
            call MPI_ALLREDUCE(&
                MPI_IN_PLACE, &
                invariant_matrix_spectrum_new(1,1,1), &
                hami%numberOfBasis*hami%numberOfBasis*hami%numberOfMatrix, &
                MPI_DOUBLE_PRECISION, &
                MPI_SUM, &
                MPI_COMM_WORLD, ierr )
    end if
    call stopCount

    !update
    deallocate(invariant_matrix_spectrum)
    allocate( invariant_matrix_spectrum &
        ( hami%numberOfBasis, hami%numberOfBasis, hami%numberOfMatrix) )
    invariant_matrix_spectrum(:,:,:) &
        = invariant_matrix_spectrum_new(:,:,:)
    deallocate(invariant_matrix_spectrum_new)

    call stopCount
end subroutine invariantMatrixForSpectrumDGEMM
