subroutine spikeZeroTemperature( iteration )
  use generall
  use parallel
  use hamiltonian, only: hami ! in
  use hamiltonian, only: basis_output ! in
  use hamiltonian, only: subspaceInfo ! in
  use hamiltonian, only: eigenvalue ! in
  use system, only: scale ! function
  use spectrum, only: positiveSpectrumOperation ! in
  use spectrum, only: negativeSpectrumOperation ! in
  use spectrum, only: binnedSpikePositive ! inout
  use spectrum, only: binnedSpikeNegative ! inout
  use spectrum, only: spectrumElementDGEMM ! function
  use spectrum, only: findSubspace ! function
  use spectrum, only: binningSpike ! function
  implicit none

  integer, intent(in) :: iteration

  integer :: i,pos, ip, iexs, iv
  integer :: countGroundState
  integer, allocatable :: groundStatePosition(:)
  integer, allocatable :: groundStateVariation(:)
  integer, allocatable :: excitedStateVariation(:)
  !the index to specify the degenerated state in ground state and excited state
  integer ::gsVariation, exsVariation
  !the index of subspace for excited state and ground state
  integer :: isubExcited, isubGround
  !GS: groundState, EXS: excited state
  !starting position of eigenvector, basis_input, number of degeneracy of subspace
  integer :: startEigenvecGS, startBasisGS, maxVariationGS
  integer ::startEigenvecEXS, startBasisEXS, maxVariationEXS

  double precision:: tmp_spike_value
  !eigenvalue of excited state
  double precision:: eigenEXS
  double precision:: tmp_spike(1,1:2)

  !partition function
  double precision :: pf
  character(50) :: spikeFile

  logical:: flag_gs, flag_exc
  integer :: ionum

  integer :: dim_GS, dim_EX
  integer :: rmax_GS, rmax_EX
  integer :: basis_min_GS, basis_min_EX
  integer :: basis_max_GS, basis_max_EX
  integer :: eigenvec_min_GS, eigenvec_min_EX
  integer :: matrixType
  double precision, allocatable :: spike_subk(:,:)

  call startCount("spikeZeroTemp")

  !find ground state and contain the index of ground state
  countGroundState=0
  do i=1, hami%numberOfBasis
     if (eigenvalue(i) < thresouldOfDegeneracy) then
        countGroundState=countGroundState+1
     end if
  end do
  pf=dble(countGroundState)

  if(my_rank .eq. 0) then
     print*, "pf=", pf
  end if

  allocate( groundStatePosition(countGroundState) )
  pos=1
  do i=1, hami%numberOfBasis
     if (eigenvalue(i) < thresouldOfDegeneracy) then
        groundStatePosition(pos)=i
        pos=pos+1
     end if
  end do

  allocate( groundStateVariation(hami%numberOfVariation) )
  allocate( excitedStateVariation(hami%numberOfVariation) )

     do ip=1, hami%numberOfMatrix
        ionum=100+ip
        write(spikeFile,'("step_",i3.2,"_itr_",i3.2,".positive_type",i3.3)') 2, iteration, ip
        open(ionum, file=trim(spikeFile))
    end do

    do ip=1, hami%numberOfMatrix
        ionum=200+ip
        write(spikeFile,'("step_",i3.2,"_itr_",i3.2,".negative_type",i3.3)') 2, iteration, ip
        open(ionum, file=trim(spikeFile))
    end do

  do i=1, countGroundState
     groundStateVariation(:) = basis_output(groundStatePosition(i))%basis(:)
     !     print*, "ground state=", groundStateVariation
     gsVariation=basis_output(groundStatePosition(i))%reference

     call findSubspace(groundStateVariation, isubGround, flag_gs)
     if(my_rank .eq. 0) then
        !        print*,"ground state index is:" ,isubGround
        if (.not.flag_gs) print*, "error, ground state does not found"
     end if

     maxVariationGS=subspaceInfo(isubGround)%dimension
     startEigenvecGS=subspaceInfo(isubGround)%count_eigenvector
     startBasisGS=subspaceInfo(isubGround)%start_input

     dim_GS=1
     rmax_GS=subspaceInfo(isubGround)%dimension
     eigenvec_min_GS=subspaceInfo(isubGround)%count_eigenvector
     basis_min_GS=subspaceInfo(isubGround)%start_input
     basis_max_GS=basis_min_GS+rmax_GS-1


     !positive region
     loop_matrix: do ip=1, hami%numberOfMatrix
        excitedStateVariation(:) &
             = groundStateVariation(:) + positiveSpectrumOperation(ip,:)
        !          print*, "excited state=",excitedStateVariation

        call findSubspace(excitedStateVariation, isubExcited, flag_exc)
        if (.not.flag_exc) then
           cycle loop_matrix
        end if
        !          print*, "excited state index is:" ,isubExcited
        maxVariationEXS=subspaceInfo(isubExcited)%dimension
        startEigenvecEXS=subspaceInfo(isubExcited)%count_eigenvector
        startBasisEXS=subspaceInfo(isubExcited)%start_input

        dim_EX=subspaceInfo(isubExcited)%dimension
        rmax_EX=subspaceInfo(isubExcited)%dimension
        eigenvec_min_EX=subspaceInfo(isubExcited)%count_eigenvector
        basis_min_EX=subspaceInfo(isubExcited)%start_input
        basis_max_EX=basis_min_EX+rmax_EX-1


        if(my_rank .eq. 0) then
           print*, "startEigenvecEXS, startEigenvecGS=", startEigenvecEXS, startBasisGS
        end if

        allocate (spike_subk(dim_EX, 1))
        spike_subk=0.0d0

        call spectrumElementDGEMM &
       ( dim_EX,rmax_EX, eigenvec_min_EX, basis_min_EX, basis_max_EX, &
         dim_GS,rmax_GS, eigenvec_min_GS, basis_min_GS, basis_max_GS, &
         ip, spike_subk)

        do iexs=startBasisEXS, startBasisEXS+maxVariationEXS-1

           eigenEXS=eigenvalue(iexs)

           tmp_spike_value = spike_subk(iexs-startBasisEXS+1,1)
           !           print*, "tmp_spike_value=", tmp_spike_value
           tmp_spike(1,1)=eigenEXS*scale(iteration)
           tmp_spike(1,2)=tmp_spike_value*tmp_spike_value/pf

                ionum=100+ip
                if (tmp_spike(1,1) >1.0D-16) then
                    write(ionum,*) tmp_spike(1,1:2)
                end if

           call binningSpike(binnedSpikePositive,tmp_spike,ip)
        end do

        deallocate (spike_subk)
     end do loop_matrix

     !negative region
     loop_matrix2:do ip=1, hami%numberOfMatrix
        excitedStateVariation(:) &
             = groundStateVariation(:) + negativeSpectrumOperation(ip,:)

        call findSubspace(excitedStateVariation,isubExcited,flag_exc)
        if (.not.flag_exc) then
           cycle loop_matrix2
        end if

        maxVariationEXS=subspaceInfo(isubExcited)%dimension
        startEigenvecEXS=subspaceInfo(isubExcited)%count_eigenvector
        startBasisEXS=subspaceInfo(isubExcited)%start_input

        dim_EX=subspaceInfo(isubExcited)%dimension
        rmax_EX=subspaceInfo(isubExcited)%dimension
        eigenvec_min_EX=subspaceInfo(isubExcited)%count_eigenvector
        basis_min_EX=subspaceInfo(isubExcited)%start_input
        basis_max_EX=basis_min_EX+rmax_EX-1

        allocate (spike_subk(1, dim_EX))
        spike_subk=0.0d0

        call spectrumElementDGEMM &
       ( dim_GS,rmax_GS, eigenvec_min_GS, basis_min_GS, basis_max_GS, &
         dim_EX,rmax_EX, eigenvec_min_EX, basis_min_EX, basis_max_EX, &
         ip, spike_subk)

        do iexs=startBasisEXS, startBasisEXS+maxVariationEXS-1
           eigenEXS=eigenvalue(iexs)

           tmp_spike_value = spike_subk(1, iexs-startBasisEXS+1)
           tmp_spike(1,1)=-eigenEXS*scale(iteration)
           tmp_spike(1,2)=tmp_spike_value*tmp_spike_value/pf

                ionum=200+ip
                if (abs(tmp_spike(1,1)) >1.0D-16) then
                    write(ionum,*) tmp_spike(1,1:2)
                end if

           call binningSpike(binnedSpikeNegative,tmp_spike,ip)
        end do
        deallocate (spike_subk)
     end do loop_matrix2
  end do

  call stopCount
end subroutine spikeZeroTemperature

subroutine spikeCFS(iteration,numberOfBasisFull)
  use generall
  use parallel
  use hamiltonian, only: hami ! in
  use hamiltonian, only: basis_output ! in
  use hamiltonian, only: subspaceInfo ! in
  use system, only: fullSubspaceInfo ! in
  use hamiltonian, only: eigenvalue ! in
  use system, only: scale ! function
  use system, only: reducedDensityMatrix ! in
  use spectrum, only: positiveSpectrumOperation ! in
  use spectrum, only: negativeSpectrumOperation ! in
  use spectrum, only: binnedSpikePositive ! inout
  use spectrum, only: binnedSpikeNegative ! inout
  use spectrum, only: spectrumElementDGEMM ! function
  use spectrum, only: findStartOutput ! function
  use spectrum, only: findSubspace ! function
  use spectrum, only: findFullSubspace ! function
  use spectrum, only: binningSpike ! function
  use spectrum, only: log10Emin,meshOfBin,numberOfBin ! in
  implicit none
  include "mpif.h"

  integer, intent(in) :: iteration
  integer, intent(in) :: numberOfBasisFull !the total number of states before truncation

  integer :: i,il, ip
  integer :: ik1, ik2,ik, imat
  integer :: ifullsubk, isubk
  integer :: maxVariationInk, keptVariationInk, reducedDensityPositionk
  integer :: variationInk, startBasis_inputk, startBasis_outputk
  integer :: startEigenveck, k1Variation, k2Variation, kVariation
  integer :: isubl
  integer :: lVariation
  integer :: maxVariationInl, startEigenvecl, startBasis_inputl

  integer :: dim_l, dim_k
  integer :: rmax_l, rmax_k
  integer :: basis_min_l, basis_min_k
  integer :: basis_max_l, basis_max_k
  integer :: eigenvec_min_l, eigenvec_min_k
  integer :: matrixType

  integer, allocatable :: variationForl(:)
  integer, allocatable :: variationFork(:)
  double precision :: El, Ek1, Ediff
  double precision :: spike_k1, spike_k2
  double precision, allocatable :: total_sum_negative(:), total_sum_positive(:)
  double precision, allocatable :: total_sum_bin_negative(:), total_sum_bin_positive(:)
  double precision :: cfsSpike(1,1:2)
  double precision, allocatable :: spike_subk(:,:)
  double precision, allocatable :: tmpBinnedSpikePositive(:,:), tmpBinnedSpikeNegative(:,:)
  !for openMP threading
  double precision, allocatable :: tmpBinnedSpikePositive_thread(:,:)
  double precision, allocatable :: tmpBinnedSpikeNegative_thread(:,:)
  double precision, allocatable :: total_sum_negative_thread(:), total_sum_positive_thread(:)


  integer :: work, work2, loadmin, loadmax

  character(50) :: spikefile
  integer :: fileIndex
  logical :: flag_l, flag_k, flag_full_k

  integer :: ierr


  call startCount("spikeCFS")
  call startCount("spikeCFS:alloc")
  if( my_rank .eq. 0 ) print*,"spike CFS start"

  !sort the basis_output


  !this subroutine should call after truncation procedure and before calculate invariant_matrix_spectrum
  !In this subroutine, I assume that
  !basis_output from 1st to numberOfBasis th element is sorted by subspace

  !loadbalancing for parallelization
  work=(numberOfBasisFull-hami%numberOfBasis)/numberOfProcess
  work2=mod(numberOfBasisFull-hami%numberOfBasis,numberOfProcess)

  !count of basis start from 1, but count of rank start from 0
  loadmin=hami%numberOfBasis+1+my_rank*work+min(my_rank,work2)
  loadmax=loadmin+work-1
  if(work2>my_rank) loadmax=loadmax+1
  !prepare temporary container for binned spike
  allocate( tmpBinnedSpikePositive(numberOfBin*hami%numberOfMatrix,2) )
  allocate( tmpBinnedSpikeNegative(numberOfBin*hami%numberOfMatrix,2) )
  allocate( tmpBinnedSpikePositive_thread(numberOfBin*hami%numberOfMatrix,2) )
  allocate( tmpBinnedSpikeNegative_thread(numberOfBin*hami%numberOfMatrix,2) )

  tmpBinnedSpikePositive=0.0d0
  tmpBinnedSpikeNegative=0.0d0
  tmpBinnedSpikePositive_thread=0.0d0
  tmpBinnedSpikeNegative_thread=0.0d0

  do imat=1, hami%numberOfMatrix
     do i=1, numberOfBin
        tmpBinnedSpikePositive(i+(imat-1)*numberOfBin,1) &
             = +10.0d0**(log10Emin+dble(i-1)/dble(meshOfBin))
        tmpBinnedSpikeNegative(i+(imat-1)*numberOfBin,1) &
             = -10.0d0**(log10Emin+dble(i-1)/dble(meshOfBin))
        tmpBinnedSpikePositive_thread(i+(imat-1)*numberOfBin,1) &
             = +10.0d0**(log10Emin+dble(i-1)/dble(meshOfBin))
        tmpBinnedSpikeNegative_thread(i+(imat-1)*numberOfBin,1) &
             = -10.0d0**(log10Emin+dble(i-1)/dble(meshOfBin))
     end do
  end do

  !loop for truncated state
  allocate( variationForl(hami%numberOfVariation) )
  allocate( variationFork(hami%numberOfVariation) )

  !initialize total sum for checking the binning procedure
  allocate( total_sum_negative(hami%numberOfMatrix) )
  allocate( total_sum_positive(hami%numberOfMatrix) )
  allocate( total_sum_negative_thread(hami%numberOfMatrix) )
  allocate( total_sum_positive_thread(hami%numberOfMatrix) )
  total_sum_negative=0.0d0
  total_sum_positive=0.0d0
  total_sum_negative_thread=0.0d0
  total_sum_positive_thread=0.0d0

  call stopCount
  call startCount("spikeCFS:omp")



  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP FIRSTPRIVATE(tmpBinnedSpikePositive_thread,tmpBinnedSpikeNegative_thread)&
  !$OMP FIRSTPRIVATE(total_sum_negative_thread,  total_sum_positive_thread)&
  !$OMP PRIVATE(i, imat)
  !$OMP DO &
  !$OMP SCHEDULE(dynamic) &
  !$OMP PRIVATE(il, ip, ik1, ik2, ik, ifullsubk, isubk, maxVariationInk, keptVariationInk)&
  !$OMP PRIVATE(reducedDensityPositionk, variationInk, startBasis_inputk, startBasis_outputk, startEigenveck)&
  !$OMP PRIVATE(k1Variation, k2Variation, kVariation, isubl, lVariation, maxVariationInl, startEigenvecl, startBasis_inputl)&
  !$OMP PRIVATE(variationForl, variationFork, El, Ek1, Ediff, spike_k1, spike_k2)&
  !$OMP PRIVATE(cfsSpike,spike_subk,flag_l,flag_k, flag_full_k)&
  !$OMP PRIVATE( dim_l, dim_k)&
  !$OMP PRIVATE( rmax_l, rmax_k)&
  !$OMP PRIVATE( basis_min_l, basis_min_k)&
  !$OMP PRIVATE( basis_max_l, basis_max_k)&
  !$OMP PRIVATE( eigenvec_min_l, eigenvec_min_k)&
  !$OMP PRIVATE( matrixType)
  do il=loadmin, loadmax
     lVariation=basis_output(il)%reference
     El=eigenvalue(il)

     variationForl(:) = basis_output(il)%basis(:)
     call findSubspace(variationForl,isubl,flag_l)

     maxVariationInl=subspaceInfo(isubl)%dimension

     dim_l=1
     rmax_l=subspaceInfo(isubl)%dimension

     startEigenvecl=subspaceInfo(isubl)%count_eigenvector
     startBasis_inputl=subspaceInfo(isubl)%start_input
     
     !shift the starting point of eigenvector block. (l is truncated basis and not the ground state of the corresponding subspace) 
     eigenvec_min_l=subspaceInfo(isubl)%count_eigenvector+rmax_l*(lVariation-1)
     basis_min_l=subspaceInfo(isubl)%start_input
     basis_max_l=basis_min_l+rmax_l-1

     print*, "iteration=", iteration, ":Gii part"

     !G_ii part
     !<k1|d^\dagger|l>*<k2|d^\dagger|l>*rhored(k2,k1)*delta(w-(Ek1-El))
     !k1 and k2 is retained state (1~numberOfBasis)
     !El > max(Ek1, Ek2) thus, w is negative
     !for determine the subspace for k1,k2, positiveSpectrumOperation can be used

     !loop for matrix variations in positiveSpectrumOperation
     loop_matrix: do ip=1, hami%numberOfMatrix
        variationFork(:) = variationForl(:) + positiveSpectrumOperation(ip,:)

        !find the reference position in "full"subspace information to determine the position of correspond reduced density matrix
        call findFullSubspace(variationFork, ifullsubk, iteration, flag_full_k)
        if (.not.flag_full_k) then
           cycle loop_matrix
        end if
        call findSubspace(variationFork,isubk,flag_k)

        startBasis_outputk=findStartOutput(variationFork )

        maxVariationInk=fullSubspaceInfo(ifullsubk)%dimension_before
        keptVariationInk=fullSubspaceInfo(ifullsubk)%dimension_after
        reducedDensityPositionk=fullSubspaceInfo(ifullsubk)%start_matrix

        variationInk=subspaceInfo(isubk)%dimension
        startEigenveck=subspaceInfo(isubk)%count_eigenvector
        startBasis_inputk=subspaceInfo(isubk)%start_input

        dim_k=fullSubspaceInfo(ifullsubk)%dimension_after
        rmax_k=subspaceInfo(isubk)%dimension
        eigenvec_min_k=subspaceInfo(isubk)%count_eigenvector
        basis_min_k=subspaceInfo(isubk)%start_input
        basis_max_k=basis_min_k+rmax_k-1

        if(maxVariationInk .ne. variationInk) then
           if(my_rank .eq. 0) then
              print*, "Warning: variation in fullSubspaceInfo and subspaceInfo is different"
           end if
           stop
        end if

        allocate( spike_subk(dim_k,1) )
        spike_subk=0.0d0
        !calculate the spike_subk using dgemm

        call spectrumElementDGEMM &
       ( dim_k,rmax_k, eigenvec_min_k, basis_min_k, basis_max_k, &
         dim_l,rmax_l, eigenvec_min_l, basis_min_l, basis_max_l, &
         ip, spike_subk)
         !print*, "subspace index for k:", isubk
         !print*, "spike_subk"
         !print*, spike_subk


        do ik1=startBasis_outputk, startBasis_outputk+keptVariationInk-1
           k1Variation=basis_output(ik1)%reference
           Ek1=eigenvalue(ik1)
           cfsSpike=0.0d0
           cfsSpike(1,1)=(Ek1-El)*scale(iteration)

           spike_k1=spike_subk(ik1-startBasis_outputk+1,1)

           do ik2=startBasis_outputk,startBasis_outputk+keptVariationInk-1
              k2Variation=basis_output(ik2)%reference
              spike_k2=spike_subk(ik2-startBasis_outputk+1,1)

              cfsSpike(1,2) = cfsSpike(1,2) + spike_k1*spike_k2 &
                   * reducedDensityMatrix &
                   ( reducedDensityPositionk+keptVariationInk &
                   *(k2Variation-1)+(k1Variation-1) )
           end do

           total_sum_negative_thread(ip) = total_sum_negative_thread(ip) + cfsSpike(1,2)

           call binningSpike(tmpBinnedSpikeNegative_thread,cfsSpike,ip)
        end do
        deallocate(spike_subk)

     end do loop_matrix

     !G_iii part
     !<l|d^\dagger|k1>*<l|d^\dagger|k2>*rhored(k2,k1)*delta(w-(El-Ek1))
     !k1 and k2 is retained state (1~numberOfBasis)
     !El > max(Ek1, Ek2) thus, w is positive --> bin to the binnedSpikePositive
     !for determine the subspace for k1,k2, negativeSpectrumOperation can be used

     !loop for matrix variations in positiveSpectrumOperation
     loop_matrix2: do ip=1, hami%numberOfMatrix
        variationFork(:) = variationForl(:) + negativeSpectrumOperation(ip,:)

        !find the reference position in "full"subspace information to determine the position of correspond reduced density matrix
        call findFullSubspace(variationFork, ifullsubk, iteration, flag_full_k)
        call findSubspace(variationFork,isubk,flag_k)

        if (.not.flag_full_k .or. .not.flag_k ) then
           cycle loop_matrix2
        end if

        startBasis_outputk=findStartOutput(variationFork )


        maxVariationInk=fullSubspaceInfo(ifullsubk)%dimension_before
        keptVariationInk=fullSubspaceInfo(ifullsubk)%dimension_after
        reducedDensityPositionk=fullSubspaceInfo(ifullsubk)%start_matrix

        variationInk=subspaceInfo(isubk)%dimension
        startEigenveck=subspaceInfo(isubk)%count_eigenvector
        startBasis_inputk=subspaceInfo(isubk)%start_input

        dim_k=fullSubspaceInfo(ifullsubk)%dimension_after
        rmax_k=subspaceInfo(isubk)%dimension
        eigenvec_min_k=subspaceInfo(isubk)%count_eigenvector
        basis_min_k=subspaceInfo(isubk)%start_input
        basis_max_k=basis_min_k+rmax_k-1

        if(maxVariationInk .ne. variationInk) then
           if(my_rank .eq. 0) then
              print*, "Warning: variation in fullSubspaceInfo and subspaceInfo is different"
           end if
           stop
        end if

        allocate( spike_subk(1,dim_k) )
        spike_subk=0.0d0
        call spectrumElementDGEMM &
       ( dim_l,rmax_l, eigenvec_min_l, basis_min_l, basis_max_l, &
         dim_k,rmax_k, eigenvec_min_k, basis_min_k, basis_max_k, &
         ip, spike_subk)


        do ik1=startBasis_outputk, startBasis_outputk+keptVariationInk-1
           k1Variation=basis_output(ik1)%reference
           Ek1=eigenvalue(ik1)

           cfsSpike=0.0d0
           cfsSpike(1,1)=(El-Ek1)*scale(iteration)
           spike_k1=spike_subk(1,ik1-startBasis_outputk+1)

           do ik2=startBasis_outputk,startBasis_outputk+keptVariationInk-1
              k2Variation=basis_output(ik2)%reference
              spike_k2=spike_subk(1,ik2-startBasis_outputk+1)
              cfsSpike(1,2) = cfsSpike(1,2) + spike_k1*spike_k2 &
                   * reducedDensityMatrix &
                   (   reducedDensityPositionk+keptVariationInk*(k2Variation-1)+(k1Variation-1))

           end do

           total_sum_positive_thread(ip) = total_sum_positive_thread(ip) + cfsSpike(1,2)
           call binningSpike(tmpBinnedSpikePositive_thread,cfsSpike,ip)
        end do

        deallocate(spike_subk)
     end do loop_matrix2
  end do
  !$OMP END DO
  !$OMP BARRIER
  !$OMP CRITICAL
  total_sum_negative=total_sum_negative+total_sum_negative_thread
  total_sum_positive=total_sum_positive+total_sum_positive_thread
  do imat=1, hami%numberOfMatrix
     do i=1, numberOfBin
        tmpBinnedSpikePositive(i+(imat-1)*numberOfBin,2) &
             = tmpBinnedSpikePositive(i+(imat-1)*numberOfBin,2) &
             + tmpBinnedSpikePositive_thread(i+(imat-1)*numberOfBin,2)
        tmpBinnedSpikeNegative(i+(imat-1)*numberOfBin,2) &
             = tmpBinnedSpikeNegative(i+(imat-1)*numberOfBin,2) &
             + tmpBinnedSpikeNegative_thread(i+(imat-1)*numberOfBin,2)
     end do
  end do
  !$OMP END CRITICAL
  !$OMP END PARALLEL

  call stopCount ! spikeCFS:omp
  call startCount("spikeCFS:sum")

  allocate( total_sum_bin_negative(hami%numberOfMatrix) )
  allocate( total_sum_bin_positive(hami%numberOfMatrix) )
  total_sum_bin_negative=0.0d0
  total_sum_bin_positive=0.0d0
  !checking the difference of sum in binned spike and native data
  !negative
  do ip=1, hami%numberOfMatrix
     do i=1, numberOfBin
        total_sum_bin_negative(ip) &
             = total_sum_bin_negative(ip) &
             + tmpBinnedSpikeNegative(i+(ip-1)*numberOfBin,2)
     end do

     if (abs(total_sum_bin_negative(ip)-total_sum_negative(ip)) &
          .ge. thresouldOfDegeneracy) then
        if(my_rank .eq. 0) then
           print*, "Warning!, sum of binned spike and native data in negative spectrum of type:", ip, &
                "is different by:", total_sum_bin_negative(ip)-total_sum_negative(ip)
        end if
     end if
  end do

  !positive
  do ip=1, hami%numberOfMatrix
     do i=1, numberOfBin
        total_sum_bin_positive(ip) &
             = total_sum_bin_positive(ip) &
             + tmpBinnedSpikePositive(i+(ip-1)*numberOfBin,2)
     end do

     if (abs(total_sum_bin_positive(ip)-total_sum_positive(ip)) &
          .ge. thresouldOfDegeneracy) then
        if(my_rank .eq. 0) then
           print*, "Warning!, sum of binned spike and native data in positive spectrum is different by:", &
                total_sum_bin_positive(ip)-total_sum_positive(ip)
        end if
     end if
  end do

  call stopCount

  !Gather and Scatter the summation by MPI_ALLREDUCE
  !MPI_ALLREDUCE seems to big rounding error??

  call startCount("spikeCFS:Allreduce")
  call MPI_ALLREDUCE( MPI_IN_PLACE, tmpBinnedSpikeNegative(1,2), &
       numberOfBin*hami%numberOfMatrix, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, tmpBinnedSpikePositive(1,2), &
       numberOfBin*hami%numberOfMatrix, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, total_sum_negative, &
       hami%numberOfMatrix, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, total_sum_positive, &
       hami%numberOfMatrix, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, total_sum_bin_negative, &
       hami%numberOfMatrix, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, total_sum_bin_positive, &
       hami%numberOfMatrix, &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
  call stopCount

  !TODO marge
  call stopCount
  call startCount("spikeCFS:accum")
  do i=1, numberOfBin*hami%numberOfMatrix
     if (binnedSpikeNegative(i,1) .ne. tmpBinnedSpikeNegative(i,1)) then
        if(my_rank .eq. 0) then
           print*, "error in merge binning, energy of tmp and total in negative region is different"
        end if
     end if
     binnedSpikeNegative(i,2) &
          = binnedSpikeNegative(i,2)+tmpBinnedSpikeNegative(i,2)

     if (binnedSpikePositive(i,1) .ne. tmpBinnedSpikePositive(i,1)) then
        if(my_rank .eq. 0) then
           print*, "error in merge binning, energy of tmp and total in negative region is different"
        end if
     end if
     binnedSpikePositive(i,2) &
          = binnedSpikePositive(i,2)+tmpBinnedSpikePositive(i,2)
  end do

  call stopCount
end subroutine spikeCFS
