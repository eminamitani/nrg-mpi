subroutine invariantMatrixForSpectrum
  use omp_lib
  use generall
  use parallel
  use hamiltonian, only: Basis_type ! type
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
  type(basis_type) :: conserv_left, conserv_right
  integer :: chargeStart, chargeEnd, spinStart, spinEnd
  integer, allocatable :: difference(:)
  logical, allocatable :: flags(:)
  integer :: imat
  logical :: flag_total
  integer :: countTrue
  integer, allocatable :: countMatrix(:)
  integer :: countIndex
  integer :: isub
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
  integer :: imatrix
  double precision :: coef

  logical :: flag_parallel !flag for parallel or not
  integer :: ierr
  integer :: p, m

  double precision, allocatable :: invariant_matrix_spectrum_new(:,:,:) !the renew value

  call startCount("invMatSpec")
  call startCount("invMatSpec:omp")
  !if the numberOfBasis < numberOfProcess, avoid to paralell
  if (hami%numberOfBasis < numberOfProcess*4 ) then
     if (my_rank .eq. 0) print*, "invariantMatrixForSpectrum:: too small matrix, avoid parallel process"
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

  allocate( invariant_matrix_spectrum_new &
       ( hami%numberOfBasis, hami%numberOfBasis, hami%numberOfMatrix) )

  invariant_matrix_spectrum_new = 0.0d0


  !$OMP PARALLEL DO &
  !$OMP SCHEDULE(dynamic) &
  !$OMP PRIVATE(ileft,iright) &
  !$OMP PRIVATE(conserv_left,conserv_right) &
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
  !$OMP PRIVATE(imatrix,coef)
  do ileft=myloadmin, myloadmax
     conserv_left = basis_output(ileft)
     
     loop_right:do iright=1,hami%numberOfBasis
        conserv_right = basis_output(iright)
        
        !calculate <ileft | d^dagger |iright> type  invariant matrix element for single particle excitation spectrum

        !except the physically unmeaning pattern,
        !this depends on the type of model Hamiltonian
        !total charge and total Sz conservation law
        chargeStart = hami%conservationBlock(1,1)
        chargeEnd   = hami%conservationBlock(1,2)
        spinStart   = hami%conservationBlock(2,1)
        spinEnd     = hami%conservationBlock(2,2)

        allocate( difference(hami%numberOfVariation) )
        allocate( flags(hami%numberOfMatrix) )

        difference(:) = conserv_left%basis(:)-conserv_right%basis(:)

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
           cycle loop_right
        end if

        !seeking the origin of each basis
        !first, find corresponding subspace
        do isub=1, numberOfSubspace
           ! count difference
           if( all(conserv_right%basis(:) == subspaceInfo(isub)%basis(:)) ) then
              rmax_right = subspaceInfo(isub)%dimension
              eigenvec_min_right = subspaceInfo(isub)%count_eigenvector
              basis_min_right = subspaceInfo(isub)%start_input
              basis_max_right = basis_min_right+rmax_right-1
           end if

           if( all(conserv_left%basis(:) == subspaceInfo(isub)%basis(:)) ) then
              rmax_left = subspaceInfo(isub)%dimension
              eigenvec_min_left = subspaceInfo(isub)%count_eigenvector
              basis_min_left = subspaceInfo(isub)%start_input
              basis_max_left = basis_min_left+rmax_left-1
           end if
        end do

        !summung up the contribution of eigenvector
        !loop for each subspace
        do isubl=basis_min_left, basis_max_left
           variation_left=basis_input(isubl)%variation
           operation_left=basis_input(isubl)%operation
           reference_left=basis_input(isubl)%reference

           loop_sub_r: do isubr=basis_min_right, basis_max_right
              variation_right=basis_input(isubr)%variation
              operation_right=basis_input(isubr)%operation
              reference_right=basis_input(isubr)%reference

              if (operation_right .ne. operation_left) then
                 cycle loop_sub_r
              end if

              diff_var_left=isubl-basis_min_left
              diff_var_right=isubr-basis_min_right

              do imat=1, size(countMatrix)
                 imatrix = countMatrix(imat)
                 coef = coefficient_invariant_matrix_spectrum &
                      (imatrix,operation_left)

                 invariant_matrix_spectrum_new(iright,ileft,imatrix) &
                      = invariant_matrix_spectrum_new(iright,ileft,imatrix) &
                      + coef &
                      * eigenvector(eigenvec_min_left &
                      +   rmax_left*(conserv_left%reference-1) &
                      +   diff_var_left) &
                      * eigenvector(eigenvec_min_right &
                      +   rmax_right*(conserv_right%reference-1) &
                      +   diff_var_right) &
                      * invariant_matrix_spectrum(reference_right,reference_left,imatrix)
              end do
           end do loop_sub_r
        end do
        deallocate(countMatrix)
     end do loop_right
  end do
  !$OMP END PARALLEL DO

  call stopCount
  call startCount("invMatSpec:Allgather")

  if (flag_parallel) then
     do m=1, hami%numberOfMatrix
        call MPI_ALLGATHERV( &
             invariant_matrix_spectrum_new(1,myloadmin,m), &
             allload(my_rank)*hami%numberOfBasis, &
             MPI_DOUBLE_PRECISION, &
             invariant_matrix_spectrum_new(1,1,m), &
             allload(:)*hami%numberOfBasis, &
             (loadmin(:)-1)*hami%numberOfBasis, &
             MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr )
     end do
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
end subroutine invariantMatrixForSpectrum
