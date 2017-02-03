subroutine postprocess( iteration )
  use generall
  use parallel
  use io_param
  use hamiltonian, only: hami ! inout
  use hamiltonian, only: Basis_type ! type
  use hamiltonian, only: basis_output ! inout
  use hamiltonian, only: basis_input ! in
  use hamiltonian, only: subspaceInfo ! in
  use hamiltonian, only: invariant_matrix_spectrum ! in
  use hamiltonian, only: eigenvalue ! inout
  use hamiltonian, only: eigenvector ! in
  use hamiltonian, only: numberOfBasis_full ! out
  use hamiltonian, only: numberOfBasis_old ! out
  use hamiltonian, only: numberOfSubspace ! in
  use hamiltonian, only: allocateType, deallocateType ! routine
  use hamiltonian, only: reallocateType, copyType ! routine
  use hamiltonian, only: lowOfInput

  use system, only: param ! in
  use system, only: mapFullSubspaceInfo ! out
  use system, only: fullEigenvalue ! out
  use system, only: fullEigenvector ! out
  use system, only: fullSubspaceInfo ! out
  use system, only: columnOfFullSubspaceInfo ! out
  use system, only: pointFullBasisInput ! out
  use system, only: pointFullBasisOutput ! out
  use system, only: pointFullEigenvector ! out
  use system, only: pointFullSubspace ! out
  use system, only: sizeOfReducedDensity ! out
  use system, only: truncationThreshold ! in
  use system, only: truncationLimitScale ! in
  use system, only: fullBasisInput ! out
  use system, only: fullBasisOutput ! out
  use system, only: allocateType, deallocateType ! routine
  use system, only: reallocateType, copyType ! routine

  use thermo, only: Cn, Sn, Xn ! out
  use thermo, only: calcThermo ! routine

  use spectrum, only: appendRestore, initRestore ! routine

  implicit none
  integer, intent(in) :: iteration

  integer :: i,j,k,ip
  integer,          allocatable :: eigenindex(:)
  type(Basis_type), allocatable :: tmp_basis_output_full(:)

  double precision, allocatable :: tmpFullEigenvector(:)
  double precision :: Ecutoff
  integer :: limit,cutoff
  integer :: degeneracy_truncation
  integer :: indexFullSubspace, startOfBasisOutput
  integer(kind=8) ::limitOfEigenvector
  character(20) :: Numitr,matrixnum
  character(50) :: format

  call startCount("postprocess")
  !calculate thermodynamical values

  call calcThermo
  if (my_rank==0) then
     write(thermoXn,*) Xn
     write(thermoSn,*) Sn
     write(thermoCn,*) Cn
  end if

  call startCount("postprocess:trunc")

  numberOfBasis_full = hami%numberOfBasis

  allocate( eigenindex(hami%numberOfBasis) )

  call sortWithIndex( hami%numberOfBasis, eigenvalue, eigenindex )

  !truncation of the basis
  !the final iteration, I avoide trancation for the CFS spectrum calculation
  if (iteration .ne. hami%indexOfFirstIteration .and. &
       iteration .ne. hami%numberOfIteration .and. &
       hami%numberOfBasis > param%truncation) then
     !set the first truncation iteration to startTrancation
     if (hami%startTrancation> iteration) then
        hami%startTrancation=iteration
        hami%flag_truncation=.true.
     end if
     limit=int(param%truncation*truncationLimitScale)
     Ecutoff=eigenvalue(param%truncation)
     i=param%truncation+1
     cutoff=param%truncation

     ! begin MIZUHO-IR
     ! modified to make sure that Ecutoff really means the cutoff energy
     Ecutoff = Ecutoff+truncationThreshold
     ! end MIZUHO-IR
     do while (i .LE. hami%numberOfBasis)
        ! begin MIZUHO-IR
        ! modified to use Ecutoff alone and compare by <= not <
        if (eigenvalue(i) <= Ecutoff ) then
!!$        if (eigenvalue(i) < Ecutoff+truncationThreshold ) then
        ! end MIZUHO-IR
           cutoff = cutoff+1
        else
           exit
        end if
        i=i+1
     end do

     do while (cutoff > limit)
        Ecutoff = Ecutoff -0.001d0
        ! begin MIZUHO-IR
        ! make the running index simple
        do i=cutoff,1,-1
           if (eigenvalue(i)>Ecutoff) then
              cutoff=cutoff-1
!!$        do i=0,cutoff-1
!!$           if (eigenvalue(cutoff-i)>Ecutoff) then
!!$              cutoff=cutoff-1
           else
              exit
           end if
        end do
     end do

     hami%numberOfBasis=cutoff
  end if

  if (my_rank .eq. 0 ) then
     !output eigenvalue
     do i=1, hami%numberOfBasis
        if(mod(iteration+1,2) .eq. 1) then
           write(eigenodd,*)  iteration+1, real(eigenvalue(i))
        else
           write(eigeneven,*) iteration+1, real(eigenvalue(i))
        end if
     end do
  end if

  !! sort in index order for keeped state after truncation
  call sortByIndex( hami%numberOfBasis, eigenvalue, eigenindex )

  !! update basis_output
  if( hami%numberOfBasis < numberOfBasis_full ) then
     call allocateType( tmp_basis_output_full, &
          numberOfBasis_full, hami%numberOfVariation )

     do i=1, numberOfBasis_full
        call copyType( tmp_basis_output_full(i), &
             basis_output(eigenindex(i)) )
     end do

     call allocateType( basis_output, &
          numberOfBasis_full, hami%numberOfVariation )

     do i=1, numberOfBasis_full
        call copyType( basis_output(i), tmp_basis_output_full(i) )
     end do

     call deallocateType( tmp_basis_output_full )
  end if

  ! out put state information of keeping state
  if (my_rank .eq. 0) then
     do i=1, hami%numberOfBasis
        write(stateInfomation,*) iteration+1, real(eigenvalue(i)), basis_output(i)%basis(:)
     end do
  end if
  deallocate(eigenindex)


  !debug
!  if (my_rank .eq. 0) then
!        !for comparing with original code
!        write(Numitr,*) iteration
!        open(110,file="itr_"//TRIM(ADJUSTL(Numitr))//"_basis_output.txt")
!
!            do i=1,numberOfBasis_full
!                write(110,*) basis_output(i)%basis(:), basis_output(i)%reference
!            end do
!
!        close(110)
!
!
!    end if

  call stopCount
  call startCount("postprocess:invariant")

  if(iteration < hami%numberOfIteration) then
     !call invariantMatrix(iteration)

     !use DGEMM version
     print*, "invariant Matrix routine"
     call invariantMatrixDGEMM(iteration)


     if (hami%flag_spectrum) call invariantMatrixForSpectrumDGEMM(iteration)

     if (my_rank .eq. 0) then
        print*,"iteration=", iteration, "size of invariant Matrix Spectrum=",size(invariant_matrix_spectrum)
     end if


  end if

  !! record the trucated number as an old one
  numberOfBasis_old=hami%numberOfBasis

  call stopCount
  call startCount("postprocess:spectrum")

  if (hami%flag_spectrum) then
     if (iteration==hami%indexOfFirstIteration) then
        call initRestore
     else
        call appendRestore( iteration )
     end if

     !in the first NRG run in cfs-NEG
     !stock the information about the basis keeping after truncation
     if(iteration .ge. hami%startTrancation) then
        if (my_rank .eq. 0) then
           print*, "truncation is done & entering the stock part"
        end if
        !containing the
        !calculate the required size of eigenvector stock

        if (iteration .eq. hami%startTrancation) then
           columnOfFullSubspaceInfo=hami%numberOfVariation+6
           pointFullBasisInput=1
           pointFullBasisOutput=1
           pointFullEigenvector=1
           pointFullSubspace=1
           allocate( mapFullSubspaceInfo(-1:hami%numberOfIteration) )
           mapFullSubspaceInfo=0
           mapFullSubspaceInfo(iteration)=1

           call allocateType &
                ( fullSubspaceInfo, numberOfSubspace, hami%numberOfVariation )
           call allocateType &
                ( fullBasisInput, numberOfBasis_full, hami%numberOfVariation )
           call allocateType &
                ( fullBasisOutput, hami%numberOfBasis, hami%numberOfVariation)
           allocate( fullEigenvalue(hami%numberOfBasis) )
           !actual size not yet clear, but here tentatively allocate the array for eigenvector
           allocate( fullEigenvector(hami%numberOfBasis*numberOfBasis_full) )

           limitOfEigenvector = hami%numberOfBasis*numberOfBasis_full
           sizeOfReducedDensity = 0
        else
           !if this is not the first truncation iteration, resize each array
           if (my_rank .eq. 0) then
              print*, "pointFullSubspace", pointFullSubspace
              print*, "numberOfSubspace", numberOfSubspace
              print*, "present size of fullSubspceInfo=", size(fullSubspaceInfo(:))
           end if

           if (pointFullSubspace+numberOfSubspace-1 > size(fullSubspaceInfo(:)) )then
              call reallocateType( fullSubspaceInfo, &
                   pointFullSubspace+numberOfSubspace-1, hami%numberOfVariation )
           end if

           call reallocateType( fullBasisInput, &
                pointFullBasisInput+numberOfBasis_full-1, hami%numberOfVariation)
           call reallocateType( fullBasisOutput, &
                pointFullBasisOutput+hami%numberOfBasis-1, hami%numberOfVariation)
           call reallocate( fullEigenvalue, &
                pointFullBasisOutput+hami%numberOfBasis-1)

           limitOfEigenvector=pointFullEigenvector+hami%numberOfBasis*numberOfBasis_full
           if (my_rank .eq. 0) then
              print*, "pointFullEigenvector::", pointFullEigenvector
              print*, "numberOfBasis", hami%numberOfBasis
              print*, "numberOfBasis_before", numberOfBasis_full
              print*, "limitOfEigenvector::",limitOfEigenvector
           end if

           !  call reallocate(fullEigenvector,pointFullEigenvector+numberOfBasis*numberOfBasis_full-1)
           !  reallocate real1D array algoritm does not match in this part

           allocate( tmpFullEigenvector(pointFullEigenvector) )
           tmpFullEigenvector(1:pointFullEigenvector) &
                = fullEigenvector(1:pointFullEigenvector)
           deallocate(fullEigenvector)
           allocate( fullEigenvector(limitOfEigenvector) )
           fullEigenvector(1:pointFullEigenvector) &
                = tmpFullEigenvector(1:pointFullEigenvector)
           deallocate(tmpFullEigenvector)

           mapFullSubspaceInfo(iteration) = pointFullSubspace
        end if
        !fillup the simple procedure for fullBasisInput, fullBasisOutput, fullEigenvalue

        do i=1, hami%numberOfBasis
           call copyType( fullBasisOutput(pointFullBasisOutput+i-1), basis_output(i) )
           fullEigenvalue(pointFullBasisOutput+i-1) = eigenvalue(i)
        end do

        do i=1, numberOfBasis_full
           call copyType( fullBasisInput(pointFullBasisInput+i-1), basis_input(i) )
        end do

        loop_subspace: do i=1, numberOfSubspace
           !fullSubspaceInfo
           !numberOfVariation+1 : the degeneracy of this subspace full truncation
           !numberOfVariation+2 : degeneracy of this subspace after truncation
           !numberOfVariation+3 : the starting point of this subspace in fullEigenvector
           !numberOfVariation+4 : the starting point of this subspace in fullBasisInput
           !numberOfVariation+5 : the starting point of this subspace in fullBasisOutput (& fullEigenvalue)
           !numberOfVariation+6 : the starting point of this subspace in reducedDensityMatrix, filled in the calculation of reduced density matrix part

           !in this loop, I fill stocks except the fullEigenvector.

           !finding the starting point in this subspace in truncated subspace
           startOfBasisOutput=hami%numberOfBasis+1
           !count the degeneracy in the keeped state in truncation
           degeneracy_truncation=0

           do j=1, hami%numberOfBasis
              if( all(subspaceInfo(i)%basis(:)==basis_output(j)%basis(:)) ) then
                 if(startOfBasisOutput > j) then
                    startOfBasisOutput=j
                 end if
                 degeneracy_truncation=degeneracy_truncation+1
              end if
           end do

           if (degeneracy_truncation > 0) then
              fullSubspaceInfo(pointFullSubspace)%basis(:) &
                   = subspaceInfo(i)%basis(:)
              fullSubspaceInfo(pointFullSubspace)%dimension_before &
                   = subspaceInfo(i)%dimension
              fullSubspaceInfo(pointFullSubspace)%dimension_after &
                   = degeneracy_truncation
              sizeOfReducedDensity &
                   = sizeOfReducedDensity &
                   + degeneracy_truncation*degeneracy_truncation
              fullSubspaceInfo(pointFullSubspace)%start_eigen &
                   = pointFullEigenvector

              fullEigenvector( pointFullEigenvector &
                   : pointFullEigenvector &
                   + subspaceInfo(i)%dimension*degeneracy_truncation-1 ) &
                   = eigenvector( subspaceInfo(i)%count_eigenvector &
                   : subspaceInfo(i)%count_eigenvector &
                   + subspaceInfo(i)%dimension &
                   * degeneracy_truncation-1 )
              !update the position index for fullEigenvector
              pointFullEigenvector = pointFullEigenvector &
                   + subspaceInfo(i)%dimension &
                   * degeneracy_truncation

              if(pointFullEigenvector .ge. limitOfEigenvector) then
                 if (my_rank .eq. 0) then
                    print*, "the size of fullEigenvector exceed the limit!"
                 end if
              end if

              fullSubspaceInfo(pointFullSubspace)%start_input &
                   = pointFullBasisInput+subspaceInfo(i)%start_input-1
              fullSubspaceInfo(pointFullSubspace)%start_output &
                   = pointFullBasisOutput+startOfBasisOutput-1
              !temporary filled by 0, actually filled in the part of reduced density matrix calculation part
              fullSubspaceInfo(pointFullSubspace)%start_matrix = 0
              !update the position index for fullSubspaceInfo
              pointFullSubspace = pointFullSubspace + 1
           end if
        end do loop_subspace

        !update the position index of each container, fullBasisInput, fullBasisOutput

        pointFullBasisInput  = pointFullBasisInput  + numberOfBasis_full
        pointFullBasisOutput = pointFullBasisOutput + hami%numberOfBasis
     end if
  end if

  call stopCount
  call stopCount
end subroutine postprocess
