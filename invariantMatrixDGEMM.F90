subroutine invariantMatrixDGEMM(iteration)
    use omp_lib
    use generall
    use parallel
    use hamiltonian, only: SubspaceInfo_type ! type
    use hamiltonian, only: hami ! in
    use hamiltonian, only: invariant_matrix ! in
    use hamiltonian, only: subspaceInfo ! in
    use hamiltonian, only: basis_output ! in
    use hamiltonian, only: basis_input ! in
    use hamiltonian, only: numberOfSubspace ! in
    use hamiltonian, only: coefficient_invariant_matrix_type ! in
    use hamiltonian, only: coefficient_invariant_matrix ! in
    use hamiltonian, only: eigenvector ! in

    implicit none
  include 'mpif.h'

  integer :: loadmin(0:numberOfProcess-1), myloadmin
  integer :: loadmax(0:numberOfProcess-1), myloadmax
  integer :: allload(0:numberOfProcess-1), myload
    integer :: iteration
    integer :: ileft, iright
    type(SubspaceInfo_type) :: subspace_left, subspace_right
    integer :: chargeStart, chargeEnd, spinStart, spinEnd
    integer :: diffQ, diffSSz
    integer :: isub, j, ip, i, ipair
    integer :: rmax_left, rmax_right
    integer :: eigenvec_min_left, eigenvec_min_right
    integer :: basis_min_left, basis_min_right
    integer :: basis_max_left, basis_max_right
    integer :: isubl, isubr
    integer :: variation_left, variation_right
    integer :: operation_left, operation_right
    integer :: diff_var_left, diff_var_right
    double precision :: coef
    integer :: matrixkind
    logical :: flag_coef
    integer :: icoef
    integer ::startl, startr
    integer, allocatable :: subspacePair(:,:)

    integer, allocatable :: keeped_basis_number(:)
    integer :: degeneracy_truncation
    !intermediate matrix for dgemm
    double precision, allocatable :: tmp1(:,:)
    !invariant matrix for subspace
    double precision, allocatable :: invariant_matrix_sub(:,:)


    double precision, allocatable :: vectorL(:,:), vectorR(:,:), coefMatrix(:,:,:)
    external :: dgemm
    double precision, parameter :: alpha=1.0, beta=0.0
    integer:: rowA, columnA, rowB, columnB, rowC, columnC

    integer :: p, m, startOfBasisOutput
    integer :: ierr
    character(50)::Numitr, matrixnum

    logical :: flag_parallel !flag for parallel or not
    integer :: IC

    call startCount("invMat")
    call startCount("invMat:omp")


    if( allocated(invariant_matrix) ) deallocate(invariant_matrix)
    allocate( invariant_matrix &
        ( hami%numberOfBasis, hami%numberOfBasis, &
        hami%numberOfConductionMatrix) )
    print*, "size of invariant_matrix:", hami%numberOfBasis

    invariant_matrix=0.0d0

    !total charge and total Sz conservation law
    chargeStart = hami%conservationBlock(1,1)
    chargeEnd   = hami%conservationBlock(1,2)
    spinStart   = hami%conservationBlock(2,1)
    spinEnd     = hami%conservationBlock(2,2)

    !construct the information of the dimension after trancation
    print*,"keeped basis"

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


    !$OMP PARALLEL DO &
    !$OMP PRIVATE(ipair, ileft,iright) &
    !$OMP PRIVATE(subspace_left,subspace_right) &
    !$OMP SHARED(chargeStart,chargeEnd,spinStart,spinEnd,keeped_basis_number) &
    !$OMP PRIVATE(diffQ,diffSSz) &
    !$OMP PRIVATE(isub) &
    !$OMP PRIVATE(rmax_left,rmax_right) &
    !$OMP PRIVATE(eigenvec_min_left,eigenvec_min_right) &
    !$OMP PRIVATE(basis_min_left,basis_min_right) &
    !$OMP PRIVATE(basis_max_left,basis_max_right) &
    !$OMP PRIVATE(isubl,isubr) &
    !$OMP PRIVATE(variation_left,variation_right) &
    !$OMP PRIVATE(operation_left,operation_right) &
    !$OMP PRIVATE(diff_var_left,diff_var_right,coef) &
    !$OMP PRIVATE(matrixkind,flag_coef,icoef)&
    !$OMP PRIVATE(vectorL, vectorR, coefMatrix, tmp1, invariant_matrix_sub)&
    !$OMP PRIVATE(startl,startr)&
    !$OMP PRIVATE(rowA, columnA, rowB, columnB, rowC, columnC)



    loop_pair:do ipair=myloadmin, myloadmax
    !$ IC=omp_get_thread_num()
        !print*, "load::", ipair, "thread::", IC

        ileft=subspacePair(ipair,1)
        iright=subspacePair(ipair,2)
        if ((keeped_basis_number(ileft) .eq. 0) .or. (keeped_basis_number(iright) .eq. 0) ) cycle loop_pair

        subspace_left = subspaceInfo(ileft)
        subspace_right = subspaceInfo(iright)

            !calculate <ileft | f^dagger |iright> type  invariant matrix element

            !except the physically unmeaning pattern,
            !this depends on the type of model Hamiltonian

            if ((chargeStart .eq. chargeEnd) .and. (chargeStart .eq. 0)) then

                !only Sz conservation

                diffSSz = sum(subspace_left%basis(spinStart:spinEnd) &
                    - subspace_right%basis(spinStart:spinEnd))
                if (abs(diffSSz) .ne. 1) then
                    cycle loop_pair
                end if

            else

                diffQ   = sum(subspace_left%basis(chargeStart:chargeEnd) &
                    - subspace_right%basis(chargeStart:chargeEnd))
                diffSSz = sum(subspace_left%basis(spinStart:spinEnd) &
                    - subspace_right%basis(spinStart:spinEnd))


                if(diffQ .ne. 1 .or. abs(diffSSz) .ne. 1) then
                    cycle loop_pair
                end if

            end if


            rmax_right         = subspaceInfo(iright)%dimension
            eigenvec_min_right = subspaceInfo(iright)%count_eigenvector
            basis_min_right    = subspaceInfo(iright)%start_input
            basis_max_right    = basis_min_right+rmax_right-1



            rmax_left          = subspaceInfo(ileft)%dimension
            eigenvec_min_left  = subspaceInfo(ileft)%count_eigenvector
            basis_min_left     = subspaceInfo(ileft)%start_input
            basis_max_left     = basis_min_left+rmax_left-1


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

            coefMatrix=0.0


            do isubl=basis_min_left, basis_max_left
                variation_left=basis_input(isubl)%variation
                operation_left=basis_input(isubl)%operation

                loop_sub_r: do isubr=basis_min_right, basis_max_right
                    variation_right=basis_input(isubr)%variation
                    operation_right=basis_input(isubr)%operation
                    if (variation_right .ne. variation_left) then
                        cycle loop_sub_r
                    end if

                    !sarching corresponding coefficient in each pair of basis
                    do icoef=1, hami%numberOfCoefficient
                        if (operation_left .eq. coefficient_invariant_matrix_type(icoef,1) &
                            .and. operation_right .eq. coefficient_invariant_matrix_type(icoef,2)) then
                            coefMatrix(isubl-basis_min_left+1,isubr-basis_min_right+1,coefficient_invariant_matrix_type(icoef,3))=coefficient_invariant_matrix(icoef)

                        end if
                    end do
                end do loop_sub_r
            end do

           ! print*, ileft,":", iright

           ! print*, "vectorL"
           ! print*, vectorL

           ! print*, "vectorR"
           ! print*, vectorR

           ! print*, "coefMatrix"

           !do ip=1, hami%numberOfConductionMatrix
           !     print*,"matrix kind:", ip
           !     print*, coefMatrix(:,:,ip)
            
           ! end do



            !dgemm

            allocate (tmp1(keeped_basis_number(ileft),rmax_right ))


            allocate (invariant_matrix_sub(keeped_basis_number(ileft), keeped_basis_number(iright)))


            do ip=1, hami%numberOfConductionMatrix

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
                    coefMatrix(:,:,ip),rowB,beta,tmp1,rowC)
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
                        invariant_matrix(startr+isubr-1, startl+isubl-1, ip)=invariant_matrix_sub(isubl,isubr)
                    end do
                end do


            end do

            deallocate(vectorL, vectorR, coefMatrix, tmp1, invariant_matrix_sub)


    end do loop_pair
     !$OMP END PARALLEL DO

    deallocate (keeped_basis_number)

    call stopCount
    call startCount("invMat:Allreduce")

    !the matrix element is zero if the calculation of matrix element was not carried out
    !and, the block of load balance in invariant matrix is not so simple
    !here I use MPI_ALLREDUCE to simplify the operation
    !MPI part is not implemented yet
!    flag_parallel=.false.
    if (flag_parallel) then
!       do m=1, hami%numberOfConductionMatrix
!            call MPI_ALLGATHERV( &
!                invariant_matrix(1,myloadmin,m), &
!                allload(my_rank)*hami%numberOfBasis, &
!                MPI_DOUBLE_PRECISION, &
!                invariant_matrix(1,1,m), &
!                allload(:)*hami%numberOfBasis, &
!                (loadmin(:)-1)*hami%numberOfBasis, &
!                MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr )
!       end do
            call MPI_ALLREDUCE(&
                MPI_IN_PLACE, &
                invariant_matrix(1,1,1), &
                hami%numberOfBasis*hami%numberOfBasis*hami%numberOfConductionMatrix, &
                MPI_DOUBLE_PRECISION, &
                MPI_SUM, &
                MPI_COMM_WORLD, ierr )
    end if
    call stopCount ! invMat:Allgather

    call stopCount

    !for debugging purpose only
!    if (my_rank .eq. 0 ) then
!        write(Numitr,*) iteration
!
!
!        do matrixkind=1, hami%numberOfConductionMatrix
!            write(matrixnum,*) matrixkind
!            open(110,file="itr_"//TRIM(ADJUSTL(Numitr))//"_invariant_matrix"//TRIM(ADJUSTL(matrixnum))//".txt")
!
!            do ileft=1, hami%numberOfBasis
!                do iright=1, hami%numberOfBasis
!                    write(110, *) invariant_matrix(iright,ileft,matrixkind)
!                end do
!            end do
!            close(110)
!        end do
!
!    end if
end subroutine invariantMatrixDGEMM
