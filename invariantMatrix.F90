subroutine invariantMatrix(iteration)
  use omp_lib
  use generall
  use parallel
  use hamiltonian, only: Basis_type ! type
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
  type(basis_type) :: conserv_left, conserv_right
  integer :: chargeStart, chargeEnd, spinStart, spinEnd
  integer :: diffQ, diffSSz
  integer :: isub
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

  integer :: p, m
  integer :: ierr
  character(50)::Numitr, matrixnum

  logical :: flag_parallel !flag for parallel or not

  call startCount("invMat")
  call startCount("invMat:omp")

  !if the numberOfBasis < numberOfProcess, avoid to paralell
  if (hami%numberOfBasis < numberOfProcess*4 ) then
     if(my_rank .eq. 0) print*, "too small matrix, avoid parallel process"
     !do calculation for all matrix in each rank
     myloadmin = 1
     myloadmax = hami%numberOfBasis
     myload    = myloadmax - myloadmin + 1
     flag_parallel= .false.
  else
     call loadbalance2( loadmin, loadmax, allload )
     myloadmin = loadmin(my_rank)
     myloadmax = loadmax(my_rank)
     myload    = allload(my_rank)
     flag_parallel= .true.
  end if

  if( allocated(invariant_matrix) ) deallocate(invariant_matrix)
  allocate( invariant_matrix &
       ( hami%numberOfBasis, hami%numberOfBasis, &
       hami%numberOfConductionMatrix) )

  invariant_matrix=0.0d0

 !total charge and total Sz conservation law
  chargeStart = hami%conservationBlock(1,1)
  chargeEnd   = hami%conservationBlock(1,2)
  spinStart   = hami%conservationBlock(2,1)
  spinEnd     = hami%conservationBlock(2,2)

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(ileft,iright) &
  !$OMP PRIVATE(conserv_left,conserv_right) &
  !$OMP SHARED(chargeStart,chargeEnd,spinStart,spinEnd) &
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
  !$OMP PRIVATE(matrixkind,flag_coef,icoef)
  do ileft=myloadmin, myloadmax
     conserv_left = basis_output(ileft)

     loop_right:do iright=1,hami%numberOfBasis
        conserv_right = basis_output(iright)

        !calculate <ileft | f^dagger |iright> type  invariant matrix element

        !except the physically unmeaning pattern,
        !this depends on the type of model Hamiltonian

        if ((chargeStart .eq. chargeEnd) .and. (chargeStart .eq. 0)) then

        !only Sz conservation

            diffSSz = sum(conserv_left%basis(spinStart:spinEnd) &
                    - conserv_right%basis(spinStart:spinEnd))
            if (abs(diffSSz) .ne. 1) then
                cycle loop_right
            end if

        else

            diffQ   = sum(conserv_left%basis(chargeStart:chargeEnd) &
                 - conserv_right%basis(chargeStart:chargeEnd))
            diffSSz = sum(conserv_left%basis(spinStart:spinEnd) &
                    - conserv_right%basis(spinStart:spinEnd))


            if(diffQ .ne. 1 .or. abs(diffSSz) .ne. 1) then
               cycle loop_right
            end if

        end if

        !seeking the origin of each basis
        !first, find corresponding subspace
        do isub=1, numberOfSubspace
           ! count difference
           if( all(conserv_right%basis(:) == subspaceInfo(isub)%basis(:)) ) then
              rmax_right         = subspaceInfo(isub)%dimension
              eigenvec_min_right = subspaceInfo(isub)%count_eigenvector
              basis_min_right    = subspaceInfo(isub)%start_input
              basis_max_right    = basis_min_right+rmax_right-1
           end if

           if( all(conserv_left%basis(:) == subspaceInfo(isub)%basis(:)) ) then
              rmax_left          = subspaceInfo(isub)%dimension
              eigenvec_min_left  = subspaceInfo(isub)%count_eigenvector
              basis_min_left     = subspaceInfo(isub)%start_input
              basis_max_left     = basis_min_left+rmax_left-1
           end if
        end do

        !summung up the contribution of eigenvector
        !loop for each subspace
        do isubl=basis_min_left, basis_max_left
           variation_left=basis_input(isubl)%variation
           operation_left=basis_input(isubl)%operation

           loop_sub_r: do isubr=basis_min_right, basis_max_right
              variation_right=basis_input(isubr)%variation
              operation_right=basis_input(isubr)%operation
              if (variation_right .ne. variation_left) then
                 cycle loop_sub_r
              end if

              diff_var_left=isubl-basis_min_left
              diff_var_right=isubr-basis_min_right
              coef=0.0d0
              flag_coef=.false.
              !sarching corresponding coefficient in each pair of basis
              do icoef=1, hami%numberOfCoefficient
                 if (operation_left .eq. coefficient_invariant_matrix_type(icoef,1) &
                      .and. operation_right .eq. coefficient_invariant_matrix_type(icoef,2)) then
                    coef=coefficient_invariant_matrix(icoef)
                    matrixkind=coefficient_invariant_matrix_type(icoef,3)
                    flag_coef=.true.
                 end if
              end do

              if (.not.flag_coef) then
                 cycle loop_sub_r
              end if

              invariant_matrix(iright,ileft,matrixkind) &
                   = invariant_matrix(iright,ileft,matrixkind) &
                   + coef &
                   * eigenvector(eigenvec_min_left &
                   +   rmax_left *(conserv_left%reference-1 ) &
                   +   diff_var_left) &
                   * eigenvector(eigenvec_min_right &
                   +   rmax_right*(conserv_right%reference-1) &
                   +   diff_var_right)
           end do loop_sub_r
        end do
     end do loop_right
  end do
  !$OMP END PARALLEL DO

  call stopCount
  call startCount("invMat:Allgather")
  if (flag_parallel) then
     do m=1, hami%numberOfConductionMatrix
        call MPI_ALLGATHERV( &
             invariant_matrix(1,myloadmin,m), &
             allload(my_rank)*hami%numberOfBasis, &
             MPI_DOUBLE_PRECISION, &
             invariant_matrix(1,1,m), &
             allload(:)*hami%numberOfBasis, &
             (loadmin(:)-1)*hami%numberOfBasis, &
             MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr )
     end do
  end if
  call stopCount ! invMat:Allgather

  call stopCount

  !for debugging purpose only
  if (my_rank .eq. 0 ) then
        write(Numitr,*) iteration


        do matrixkind=1, hami%numberOfConductionMatrix
            write(matrixnum,*) matrixkind
            open(110,file="itr_"//TRIM(ADJUSTL(Numitr))//"_invariant_matrix"//TRIM(ADJUSTL(matrixnum))//".txt")

            do ileft=1, hami%numberOfBasis
            do iright=1, hami%numberOfBasis
            write(110, *) invariant_matrix(iright,ileft,matrixkind)
            end do
            end do
            close(110)
        end do

    end if
end subroutine invariantMatrix
