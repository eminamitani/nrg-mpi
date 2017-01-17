subroutine read_initial_basis
  !************************************************
  !this subroutine read basis information from file
  !basis_input.dat : initial basis_input data, all integer

  use hamiltonian
  use io_param
  use generall
  use system
  use spectrum

  implicit none
  integer :: i,j,k,ii,imat
  double precision :: tmp_value
  integer :: tmp_index, counter
  character (50) ::tmp_character

  !this is a set up for the calculation without impurity
  !for calculate the contribution of conduction chain part on Xn, Sn etc
  !the empty state is always  the first component in the input file
  if (param%noimpurity) then
     hami%numberOfBasis=1
     sizeOfEigenvector=1
     numberOfSubspace=1
  end if

  open(20, file="initial.dat")
  read (20,*) tmp_character
  !!! print*, "read ::", tmp_character
  !!! print*, "numberOfBasis:" , hami%numberOfBasis
  do i=1, hami%numberOfBasis
     !eigenvalue is scaled already in Java code
     read(20,*) &
          basis_output(i)%basis(:), &
          basis_output(i)%reference, &
          eigenvalue(i)
     basis_output(i)%variation = 0
     basis_output(i)%operation = 0
  end do

  read (20,*) tmp_character
  !!! print*, "read ::", tmp_character

  !counter to the indicate the position in invariant_matrix
  counter =1
  do i=1, hami%numberOfConductionMatrix
     do j=1, hami%numberOfBasis
        do k=1, hami%numberOfBasis
           read(20,*) invariant_matrix(k,j,i)
!!!!           print*, counter, invariant_matrix(k,j,i)
           counter=counter+1
        end do
     end do
  end do

  read (20,*) tmp_character
  !!! print*, "read ::", tmp_character
  do i=1, hami%numberOfCoefficient
     read(20,*) coefficient_invariant_matrix_type(i,1),&
          coefficient_invariant_matrix_type(i,2),&
          coefficient_invariant_matrix_type(i,3),&
          coefficient_invariant_matrix(i)
  end do

  read (20,*) tmp_character
  !!! print*, "read ::", tmp_character
  do i=1, hami%numberOfCoefficient
     read(20,*) coefficient_diagonalization_type(i,1),&
          coefficient_diagonalization_type(i,2),&
          coefficient_diagonalization_type(i,3),&
          coefficient_diagonalization_type(i,4),&
          coefficient_diagonalization(i)
  end do

  read (20,*) tmp_character
  !!! print*, "read ::", tmp_character
  do i=1, hami%numberOfOperation
     read(20,*) operation(i,1:hami%numberOfVariation+hami%numberOfChain)
  end do

  read (20,*) tmp_character
  !!! print*, "read ::", tmp_character

  do i=1, hami%numberOfChain
     read (20,*) tmp_character
     !!! print*, "read ::", tmp_character
     do j=0, hami%numberOfIteration+1
        read(20,*) chain_diagonal(j,i)
!!!!        print*, j,  chain_diagonal(j,i)
     end do
  end do

  do i=1, hami%numberOfChain
     read (20,*) tmp_character
     !!! print*, "read ::", tmp_character
     do j=hami%indexOfFirstIteration, hami%numberOfIteration+1
        read(20,*) chain_nondiagonal(j,i)
!!!!        print*, j,  chain_nondiagonal(j,i)
     end do
  end do

  !scaling the wilson chain
  do i=1, hami%numberOfChain
     select case (param%scalingType)
     case(0)
        !disc=Wilcon&Campo no scaling
        !not tested yet (20130208)
        do ii=0, hami%numberOfIteration+1
           !disc=Wilson&Campo no scaling
           chain_diagonal(ii,i)=chain_diagonal(ii,i) &
                *param%Lambda**(dble(ii-1)*0.5d0)
           chain_nondiagonal(ii,i)=chain_nondiagonal(ii,i) &
                *param%Lambda**(dble(ii)*0.5d0)
        end do

     case(2)
        !disc=Campo
        do ii=0, hami%numberOfIteration+1
           chain_diagonal(ii,i)=chain_diagonal(ii,i) &
                *param%Lambda**(dble(ii-1)*0.5d0+param%ztrick) &
                *log(param%Lambda)/param%bandWidth/(param%Lambda-1.0d0)
           chain_nondiagonal(ii,i)=chain_nondiagonal(ii,i) &
                *param%Lambda**(dble(ii)*0.5d0+param%ztrick) &
                *log(param%Lambda)/param%bandWidth/(param%Lambda-1.0d0)
        end do

     case(1)
        !disc=Wilson scaling
        do ii=0, hami%numberOfIteration+1
           chain_diagonal(ii,i)=chain_diagonal(ii,i)&
                *param%Lambda**(dble(ii-1)*0.5d0+param%ztrick-1.0d0) &
                /param%bandWidth*2.0d0/(1.0d0+1.0d0/param%Lambda)
           chain_nondiagonal(ii,i)=chain_nondiagonal(ii,i)&
                *param%Lambda**(dble(ii)*0.5d0+param%ztrick-1.0d0) &
                /param%bandWidth*2.0d0/(1.0d0+1.0d0/param%Lambda)
           print*,"iteration",ii,"scaling factor=", &
                param%Lambda**(dble(ii)*0.5d0+param%ztrick-1.0d0) &
                /param%bandWidth*2.0d0/(1.0d0+1.0d0/param%Lambda)
        end do

     end select
  end do

  if (hami%flag_spectrum)  then
     read(20,*) tmp_character
     !!! print*, "read::", tmp_character
     do imat=1, hami%numberOfMatrix
        do i=1, hami%numberOfOperation
           read(20,*) coefficient_invariant_matrix_spectrum(imat,i)
        end do
     end do
     read(20,*) tmp_character
     !!! print*, "read::", tmp_character
     counter=1

     do imat=1, hami%numberOfMatrix
        do j=1, hami%numberOfBasis
           do i=1, hami%numberOfBasis
              read(20,*) invariant_matrix_spectrum(i,j,imat)
              counter=counter+1
           end do
        end do
     end do
     read(20,*) tmp_character
     !!! print*, "read::", tmp_character
     do i=1, hami%numberOfMatrix
        read(20,*) conservation_difference_spectrum(i,1:hami%numberOfVariation)
     end do

     read(20,*) tmp_character
     !!! print*, "read::", tmp_character

     do i=1, hami%numberOfMatrix
        read(20,*) positiveSpectrumOperation(i,1:hami%numberOfVariation)
     end do

     read(20,*) tmp_character
     !!! print*, "read::", tmp_character

     do i=1, hami%numberOfMatrix
        read(20,*) negativeSpectrumOperation(i,1:hami%numberOfVariation)
     end do
  end if

  if (hami%flag_BCS) then
!chain for BCS
    do i=1, hami%numberOfChain
    read (20,*) tmp_character
    print*, "read ::", tmp_character
    do j=0, hami%numberOfIteration+1
    read(20,*) chain_bcs(j,i)
    print*, j,  chain_bcs(j,i)
    end do
    end do
!coefficient for BCS

    read (20,*) tmp_character
    print*, "read ::", tmp_character
    do i=1, hami%numberOfBCSCoefficient
    read(20,*) coefficient_BCS(i,1),coefficient_BCS(i,2),coefficient_BCS(i,3)
    end do
  end if

  close(20)

end subroutine read_initial_basis
