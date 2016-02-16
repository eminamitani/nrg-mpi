program main
  use parallel
  use system
  use io_param
  use hamiltonian
  use spectrum
  use restart

  implicit none
  include 'mpif.h'

  !input parameter namelist
  integer :: ierr, i, imat
  integer :: iteration, iteration_restart, mode_restart
  character(50) :: spikefile

  !MPI process initialize
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, numberOfProcess, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank,ierr)

  call initCount
  call startCount("main")
  call startCount("main:setup")

  if(my_rank .eq. 0) print*, "*****MPI version of NRG code******"

  !setup Hamiltonian information
  call setupHamiltonian

  if(hami%flag_spectrum) call setSpectrumInfo

  !preperation of the container
  !20130308 fix the length of chain
  allocate( chain_diagonal(0:hami%numberOfIteration+1, hami%numberOfChain) )
  allocate( chain_nondiagonal(hami%indexOfFirstIteration:hami%numberOfIteration+1, &
       hami%numberOfChain))

  !reading most initial input parameters
  if (my_rank .eq. 0) then
     call read_parameter
  end if

  !open files to event log
  if (my_rank .eq. 0) then
     open(errorlog,file=errorfile)
     open(eigeneven,file=fileEigenEven)
     open(eigenodd,file=fileEigenOdd)
     open(thermoXn, file=fileXn)
     open(thermoSn, file=fileSn)
     open(thermoCn, file=fileCn)
     open(stateInfomation, file=fileInfo )
  end if

  !read input data from file
  if (my_rank .eq. 0) then
     call read_initial_basis
  end if

  call stopCount ! main:setup
  call startCount("main:Bcast")
  !broadcast
  call MPI_BCAST(param,sizeof(param),MPI_BYTE,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_BCAST(scalingType,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_BCAST(truncation,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_BCAST(Lambda,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_BCAST(ztrick,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!$  call MPI_BCAST(bandWidth,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  !chain_diagonal & chain_nondiagonal start from 0, keep in mind
  call MPI_BCAST( chain_diagonal(0,1), &
       size(chain_diagonal), MPI_DOUBLE_PRECISION, &
       0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( chain_nondiagonal(hami%indexOfFirstIteration,1), &
       size(chain_nondiagonal), MPI_DOUBLE_PRECISION, &
       0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( locationType(basis_output), &
       sizeofType(basis_output), MPI_BYTE, &
       0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( eigenvalue, &
       size(eigenvalue), MPI_DOUBLE_PRECISION, &
       0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( coefficient_invariant_matrix(1), &
       size(coefficient_invariant_matrix), MPI_DOUBLE_PRECISION, &
       0,MPI_COMM_WORLD, ierr )
  call MPI_BCAST( coefficient_invariant_matrix_type(1,1), &
       size(coefficient_invariant_matrix_type), MPI_INTEGER, &
       0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( invariant_matrix(1,1,1), &
       size(invariant_matrix), MPI_DOUBLE_PRECISION, &
       0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( operation(1,1), &
       size(operation), MPI_INTEGER, &
       0, MPI_COMM_WORLD,ierr )
  call MPI_BCAST( coefficient_diagonalization(1), &
       size(coefficient_diagonalization), MPI_DOUBLE_PRECISION, &
       0, MPI_COMM_WORLD, ierr )
  call MPI_BCAST( coefficient_diagonalization_type(1,1), &
       size(coefficient_diagonalization_type), MPI_INTEGER, &
       0, MPI_COMM_WORLD, ierr )

  if (hami%flag_spectrum) then
     call MPI_BCAST( coefficient_invariant_matrix_spectrum(1,1), &
          size(coefficient_invariant_matrix_spectrum), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
     call MPI_BCAST( positiveSpectrumOperation(1,1), &
          size(positiveSpectrumOperation), &
          MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
     call MPI_BCAST( negativeSpectrumOperation(1,1), &
          size(negativeSpectrumOperation), &
          MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
     call MPI_BCAST( conservation_difference_spectrum(1,1), &
          size(conservation_difference_spectrum), &
          MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
     call MPI_BCAST( invariant_matrix_spectrum(1,1,1), &
          size(invariant_matrix_spectrum), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  end if
  call stopCount ! main:bcast

  !stock the initial state
  call stockInitialState
  !calculate first invariant_matrix
  !the calculation done here is the invariant matrix for impurity site,
  !so the invariantMatrixImpurity is called
  !2014 01 09 TODO: change to read from file
  !call invariantMatrixImpurity

  !prepare the container for binnedSpike
  Emax=scale(-1)
  log10Emax=int(log10(Emax))+1.0d0
  numberOfBin=int(log10Emax-log10Emin)*meshOfBin+1

  allocate( binnedSpikePositive(numberOfBin*hami%numberOfMatrix,2) )
  allocate( binnedSpikeNegative(numberOfBin*hami%numberOfMatrix,2) )

  binnedSpikePositive=0.0d0
  binnedSpikeNegative=0.0d0

  do imat=1, hami%numberOfMatrix
     do i=1, numberOfBin
        binnedSpikePositive(i+(imat-1)*numberOfBin,1) = +10.0d0**(log10Emin+dble(i-1)/dble(meshOfBin))
        binnedSpikeNegative(i+(imat-1)*numberOfBin,1) = -10.0d0**(log10Emin+dble(i-1)/dble(meshOfBin))
     end do
  end do

  !copy this initial set up for **_old container
  numberOfBasis_old=hami%numberOfBasis
!!$  !for checking the startinc iteration of truncation
!!$  hami%startTrancation=hami%numberOfIteration+1
!!$  hami%flag_truncation= .false.

  if( param%loadRestart ) then
     call loadRestart( iteration_restart, mode_restart )
  else
     mode_restart = MODE_RESTART_NON
  end if

  if( mode_restart == MODE_RESTART_CFS ) then
     ! skip NRG
  else
     call startCount("main:nrg")
     do iteration=hami%indexOfFirstIteration, hami%numberOfIteration
        call nrgIteration(iteration)
     end do
     call stopCount ! main:setup
  end if

  !CFS-NRG
  if (hami%flag_spectrum) then
     call startCount("main:spec")

     if( mode_restart == MODE_RESTART_CFS ) then
        ! skip initialize CFS
     else
        if (hami%flag_truncation) then
           call initializeReducedDensity
        end if

        !keep the size of each container
        elementNumFullBasisInput  = pointFullBasisInput-1
        elementNumFullBasisOutput = pointFullBasisOutput-1
        elementNumFullSubspace    = pointFullSubspace-1
        elementNumReducedDensity  = pointReducedDensity-1
        elementNumFullEigenVector = pointFullEigenvector-1

        !calculate reducedDensityMatrix
        call reducedDensityMatrixCalculation
     end if

     !spectrum calculation by using CFS method
     do iteration=hami%startTrancation, hami%numberOfIteration
        if( mode_restart == MODE_RESTART_CFS ) then
           if( iteration < iteration_restart ) then
              ! skip steps before the restart file
              cycle
           end if
           if( iteration > iteration_restart ) then
              ! overwrite restart file after the step of the restart file
              if( param%saveRestart ) then
                 call saveRestart( iteration )
              end if
           end if
        else
           ! overwrite restart file at any step
           if( param%saveRestart ) then
              call saveRestart( iteration )
           end if
        end if

        call cfsSpectrumCalculation( iteration )
     end do
     call stopCount ! main:spec

     call startCount("main:output")
     if (my_rank .eq. 0) then
        do imat=1,hami%numberOfMatrix
           write(spikefile,'("CFS-spectrum-positive-type",i3.3)')  imat
           open(210,file=trim(spikefile))
           do  i=1, numberOfBin
              if(abs(binnedSpikePositive(i+(imat-1)*numberOfBin,2)) > &
                   1.0D-16*scale(hami%numberOfIteration)) then
                 write(210,*) binnedSpikePositive(i+(imat-1)*numberOfBin,1:2)
              end if
           end do
           close(210)
        end do

        do imat=1,hami%numberOfMatrix
           write(spikefile,'("CFS-spectrum-negative-type",i3.3)')  imat
           open(210,file=trim(spikefile))
           do  i=1, numberOfBin
              if(abs(binnedSpikeNegative(i+(imat-1)*numberOfBin,2)) > &
                   1.0D-16*scale(hami%numberOfIteration)) then
                 write(210,*) binnedSpikeNegative(i+(imat-1)*numberOfBin,1:2)
              end if
           end do
           close(210)
        end do
     end if
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     call stopCount ! main:output
  end if

  if (my_rank .eq. 0) then
     print*, "finalize the calculation"
  end if

  call stopCount ! main
  call showCount

  call MPI_FINALIZE(ierr)
end program main
