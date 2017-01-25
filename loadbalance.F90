
subroutine loadbalance(loadmin, loadmax)
  !very symple block separation
  !set up lower limit and upper limit in each rank
  !used in calculation of invariant matrix element
  use hamiltonian
  use parallel
  implicit none
  integer, intent(inout) :: loadmin,loadmax
  integer :: work,work2

  work=hami%numberOfBasis/numberOfProcess
  work2=mod(hami%numberOfBasis,numberOfProcess)

  !count of basis start from 1, but count of rank start from 0
  loadmin=1+my_rank*work+min(my_rank,work2)
  loadmax=loadmin+work-1
  if(work2>my_rank) loadmax=loadmax+1
end subroutine loadbalance

!simplify
subroutine loadbalanceDiagonalization(loadParRank)
  use hamiltonian
  use parallel
  use io_param

  implicit none
  integer, intent(inout) ::loadParRank(0:numberOfProcess-1,2)
  integer :: i,j, loadmin, loadmax
  integer :: work, work2

  work=numberOfSubspace/numberOfProcess
  work2=mod(numberOfSubspace,numberOfProcess)

  !count of basis start from 1, but count of rank start from 0
  do i=0, numberOfProcess-1
     loadmin=1+i*work+min(i,work2)
     loadmax=loadmin+work-1
     if(work2>i) loadmax=loadmax+1
     loadParRank(i,1)=loadmin
     loadParRank(i,2)=loadmax
  end do
end subroutine loadbalanceDiagonalization

subroutine loadbalance2( loadmin, loadmax, allload )
  use hamiltonian
  use parallel
  implicit none
  integer, intent(out) :: loadmin(0:numberOfProcess-1)
  integer, intent(out) :: loadmax(0:numberOfProcess-1) 
  integer, intent(out) :: allload(0:numberOfProcess-1) 
  integer :: p

  allload(:) = hami%numberOfBasis/numberOfProcess

  do p=0, numberOfProcess-1
     if( sum(allload(:)) /= hami%numberOfBasis ) then
        allload(p) = allload(p) + 1
     end if
  end do

  do p=0, numberOfProcess-1
     if( p==0 ) then
        loadmin(p) = 1
     else
        loadmin(p) = loadmin(p-1) + allload(p-1)
     end if
     loadmax(p) = loadmin(p) + allload(p) - 1
  end do
end subroutine loadbalance2


!subroutine for evaluate loadbalance in invariantMatrix calculation
subroutine loadbalance3( loadmin, loadmax, allload, total )
  use hamiltonian
  use parallel
  implicit none
  integer, intent(in)  :: total
  integer, intent(out) :: loadmin(0:numberOfProcess-1)
  integer, intent(out) :: loadmax(0:numberOfProcess-1)
  integer, intent(out) :: allload(0:numberOfProcess-1)
  integer :: p

  allload(:) = total/numberOfProcess

  do p=0, numberOfProcess-1
     if( sum(allload(:)) /= total ) then
        allload(p) = allload(p) + 1
     end if
  end do

  do p=0, numberOfProcess-1
     if( p==0 ) then
        loadmin(p) = 1
     else
        loadmin(p) = loadmin(p-1) + allload(p-1)
     end if
     loadmax(p) = loadmin(p) + allload(p) - 1
  end do
end subroutine loadbalance3

subroutine loadbalanceDiagonalization2( ismysubspace )
  use hamiltonian
  use parallel
  use io_param

  implicit none
  logical, intent(out) :: ismysubspace(numberOfSubspace,0:numberOfProcess-1)

  integer          :: i
  double precision :: load(numberOfSubspace)
  double precision :: loadmax
  integer          :: j, jmax
  double precision :: busy(0:numberOfProcess-1)
  double precision :: busymin
  integer          :: p, pmin


  busy(:) = 0.0d0
  ismysubspace(:,:) = .false.

  !! estimate load of each job roughly and record loads to the list
  do j=1, numberOfSubspace
     load(j) = dble(subspaceInfo(j)%dimension)**3
  end do

  !! assign all jobs to any of ranks
  do i=1, numberOfSubspace
     !! find the most heavy job left in the list
     do j=1, numberOfSubspace
        if( j==1 .or. loadmax < load(j) ) then
           loadmax = load(j)
           jmax = j
        end if
     end do

     !! find the least busy process
     do p=0, numberOfProcess-1
        if( p==0 .or. busymin > busy(p) ) then
           busymin = busy(p)
           pmin = p
        end if
     end do

     !! add the load to the process
     busy(pmin) = busy(pmin) + load(jmax)
     !! remove the load from the list
     load(jmax) = 0.0d0

     !! assign this jobs to this process
     ismysubspace(jmax,pmin) = .true.
  end do

!!$  !! for check
!!$  if( my_rank == 0 ) then
!!$     do j=1, numberOfSubspace
!!$        load(j) = dble(subspaceInfo(j)%dimension)**3
!!$     end do
!!$
!!$     write(*,"('---',$)")
!!$     do p=0, numberOfProcess-1
!!$        write(*,"(i10$)") p
!!$     end do
!!$     write(*,*)
!!$
!!$     do j=1, numberOfSubspace
!!$        write(*,"(i3$)") j
!!$        do p=0, numberOfProcess-1
!!$           if( ismysubspace(j,p) ) then
!!$              write(*,"(f10.1$)") load(j)
!!$           else
!!$              write(*,"(f10.1$)") 0.0d0
!!$           end if
!!$        end do
!!$        write(*,*)
!!$     end do
!!$
!!$     write(*,"('sum',$)")
!!$     do p=0, numberOfProcess-1
!!$        write(*,"(f10.1$)") busy(p)
!!$     end do
!!$     write(*,*)
!!$     flush(6)
!!$  end if

end subroutine loadbalanceDiagonalization2


subroutine loadbalance1DSimple(loadParRank,dimload)
  use hamiltonian
  use parallel
  use io_param

  implicit none
  integer :: dimload
  integer :: i, work, work2,loadmax,loadmin
  integer,intent(inout):: loadParRank(0:numberOfProcess-1,2)

  if (dimload < numberOfProcess) then
     print*, "serious error, too many process for small dimension of the problem"
  end if

  work=dimload/numberOfProcess
  work2=mod(dimload,numberOfProcess)

  !count of basis start from 1, but count of rank start from 0
  do i=0,numberOfProcess-1
     loadmin=1+i*work+min(i,work2)
     loadmax=loadmin+work-1
     if(work2>i) loadmax=loadmax+1
     loadParRank(i,1)=loadmin
     loadParRank(i,2)=loadmax
  end do

  write(buglog,*) "loadbalancefor1D"
  do i=0, numberOfProcess-1
     write(buglog,*) loadParRank(i,1:2)
  end do
end subroutine loadbalance1Dsimple
