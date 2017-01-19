module count_module
  implicit none

  type Count_type
     character(len=64) :: name
     real(8) :: time_started
     real(8) :: time_elapsed
     integer(8) :: count_call
     integer(8) :: count_load
     integer(8) :: count_store
     integer(8) :: count_add
     integer(8) :: count_mul
     integer(8) :: count_div
     integer(8) :: count_func
  end type Count_type

  integer, parameter, private :: ncount = 128 ! maximum counters in the program
  type(Count_type), save, private :: vcount(ncount)
  integer, private :: itarget = 0
  integer, parameter, private :: ntarget = 16 ! maximum counters at the same time
  integer, private :: vtarget(ntarget)
  integer, private :: mtarget = 0

  public :: initCount
  public :: showCount
  public :: startCount
  public :: stopCount
  public :: countLoad
  public :: countStore
  public :: countAdd
  public :: countMul
  public :: countDiv
  public :: countFunc

contains
  subroutine initCount
    implicit none
    integer :: id

    do id=1, ncount
       vcount(id)%name = ""
       vcount(id)%time_started = 0.0d0
       vcount(id)%time_elapsed = 0.0d0
       vcount(id)%count_call = 0
       vcount(id)%count_load = 0
       vcount(id)%count_store = 0
       vcount(id)%count_add = 0
       vcount(id)%count_mul = 0
       vcount(id)%count_div = 0
       vcount(id)%count_func = 0
    end do
    itarget = 0
    mtarget = 0

    return
  end subroutine initCount


  subroutine showCount
    implicit none
    include "mpif.h"
    integer :: id
    integer :: mpi_info, mpi_rank, mpi_size

    real(8) :: dwork
    integer(8) :: iwork
    integer, parameter :: unit = 100
    character(len=64) :: filename

    call MPI_Comm_size( MPI_COMM_WORLD, mpi_size, mpi_info )
    call MPI_Comm_Rank( MPI_COMM_WORLD, mpi_rank, mpi_info )

    do id=1, ncount
       if( vcount(id)%name == "" ) cycle

       dwork = vcount(id)%time_elapsed
       call MPI_Reduce( dwork, vcount(id)%time_elapsed, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_COMM_WORLD, mpi_info )
       vcount(id)%time_elapsed = vcount(id)%time_elapsed/mpi_size

       iwork = vcount(id)%count_call
       call MPI_Reduce( iwork, vcount(id)%count_call, 1, MPI_INTEGER8, &
            MPI_SUM, 0, MPI_COMM_WORLD, mpi_info )
       vcount(id)%count_call = vcount(id)%count_call/mpi_size

       iwork = vcount(id)%count_load
       call MPI_Reduce( iwork, vcount(id)%count_load, 1, MPI_INTEGER8, &
            MPI_SUM, 0, MPI_COMM_WORLD, mpi_info )
       vcount(id)%count_load = vcount(id)%count_load/mpi_size

       iwork = vcount(id)%count_store
       call MPI_Reduce( iwork, vcount(id)%count_store, 1, MPI_INTEGER8, &
            MPI_SUM, 0, MPI_COMM_WORLD, mpi_info )
       vcount(id)%count_store = vcount(id)%count_store/mpi_size

       iwork = vcount(id)%count_add
       call MPI_Reduce( iwork, vcount(id)%count_add, 1, MPI_INTEGER8, &
            MPI_SUM, 0, MPI_COMM_WORLD, mpi_info )
       vcount(id)%count_add = vcount(id)%count_add/mpi_size

       iwork = vcount(id)%count_mul
       call MPI_Reduce( iwork, vcount(id)%count_mul, 1, MPI_INTEGER8, &
            MPI_SUM, 0, MPI_COMM_WORLD, mpi_info )
       vcount(id)%count_mul = vcount(id)%count_mul/mpi_size

       iwork = vcount(id)%count_div
       call MPI_Reduce( iwork, vcount(id)%count_div, 1, MPI_INTEGER8, &
            MPI_SUM, 0, MPI_COMM_WORLD, mpi_info )
       vcount(id)%count_div = vcount(id)%count_div/mpi_size

       iwork = vcount(id)%count_func
       call MPI_Reduce( iwork, vcount(id)%count_func, 1, MPI_INTEGER8, &
            MPI_SUM, 0, MPI_COMM_WORLD, mpi_info )
       vcount(id)%count_func = vcount(id)%count_func/mpi_size
    end do
    
    if( mpi_rank == 0 ) then
       write(filename,"('count.rank',i1,'.csv')") mpi_rank
       filename="count.csv"
       open(unit,file=filename)

       write(unit,'("# -----------------------------------------")')
       write(unit,'("# routine,  time[sec]/process/thread,  call")') 
!!$       write(unit,'("# file-routine-loop, time[sec],    call,    load,   store,     add,     mul,     div,    func")') 

       do id=1, ncount
          if( vcount(id)%name == "" ) cycle

          write(unit,'(a,",",f10.3,",",i8)') &
               trim(vcount(id)%name), &
               vcount(id)%time_elapsed, &
               vcount(id)%count_call

!!$          write(unit,'(a,",",f10.3,",",i8,6(",",i8))') &
!!$               trim(vcount(id)%name), &
!!$               vcount(id)%time_elapsed, &
!!$               vcount(id)%count_call, &
!!$               vcount(id)%count_load, &
!!$               vcount(id)%count_store, &
!!$               vcount(id)%count_add, &
!!$               vcount(id)%count_mul, &
!!$               vcount(id)%count_div, &
!!$               vcount(id)%count_func
       end do
       write(unit,'("# -----------------------------------------")')
       close(unit)
    end if

    return
  end subroutine showCount



  subroutine startCount( name )
    implicit none
    character(len=*) :: name
    include "mpif.h"
    integer :: id

    !$OMP MASTER
    itarget = 0
    do id=1, ncount
       if( vcount(id)%name == name ) then
          itarget = id
          exit
       end if
    end do

    if( itarget == 0 ) then
       do id=1, ncount
          if( vcount(id)%name == "" ) then
             itarget = id
             exit
          end if
       end do
    else
       if( vcount(itarget)%time_started > 0.0d0 ) then
          write(*,*) "Error: count has already started for ", &
               trim(vcount(itarget)%name)
          stop
       end if
    end if

    mtarget = mtarget + 1
    vtarget(mtarget) = itarget

    vcount(itarget)%name = name
    vcount(itarget)%time_started = MPI_wtime()
    vcount(itarget)%count_call &
         = vcount(itarget)%count_call + 1

!!$    write(*,*) "DEBUG: start ", trim(vcount(itarget)%name)
    !$OMP END MASTER
  end subroutine startCount

  subroutine stopCount
    implicit none
    include "mpif.h"

    !$OMP MASTER
    if( itarget == 0 ) then
       write(*,*) "Error: count has not started yet."
       stop
    end if

!!$    write(*,*) "DEBUG: stop  ", trim(vcount(itarget)%name)

    vcount(itarget)%time_elapsed &
         = vcount(itarget)%time_elapsed &
         + ( MPI_wtime() - vcount(itarget)%time_started )
    vcount(itarget)%time_started = 0.0d0

    mtarget = mtarget - 1
    if( mtarget > 0 ) then
       itarget = vtarget(mtarget)
    else
       itarget = 0
    end if
    !$OMP END MASTER
  end subroutine stopCount

  subroutine countLoad( count )
    implicit none
    integer, intent(in) :: count

    if( itarget == 0 ) then
       write(*,*) "Error: count has not started yet."
       stop
    end if

    !$OMP CRITICAL
    vcount(itarget)%count_load &
         = vcount(itarget)%count_load + count
    !$OMP END CRITICAL

    return
  end subroutine countLoad

  subroutine countStore( count )
    implicit none
    integer, intent(in) :: count

    if( itarget == 0 ) then
       write(*,*) "Error: count has not started yet."
       stop
    end if

    !$OMP CRITICAL
    vcount(itarget)%count_store &
         = vcount(itarget)%count_store + count
    !$OMP END CRITICAL

    return
  end subroutine countStore

  subroutine countAdd( count )
    implicit none
    integer, intent(in) :: count

    if( itarget == 0 ) then
       write(*,*) "Error: count has not started yet."
       stop
    end if

    !$OMP CRITICAL
    vcount(itarget)%count_add &
         = vcount(itarget)%count_add + count
    !$OMP END CRITICAL

    return
  end subroutine countAdd

  subroutine countMul( count )
    implicit none
    integer, intent(in) :: count

    if( itarget == 0 ) then
       write(*,*) "Error: count has not started yet."
       stop
    end if

    !$OMP CRITICAL
    vcount(itarget)%count_mul &
         = vcount(itarget)%count_mul + count
    !$OMP END CRITICAL

    return
  end subroutine countMul

  subroutine countDiv( count )
    implicit none
    integer, intent(in) :: count

    if( itarget == 0 ) then
       write(*,*) "Error: count has not started yet."
       stop
    end if

    !$OMP CRITICAL
    vcount(itarget)%count_div &
         = vcount(itarget)%count_div + count
    !$OMP END CRITICAL

    return
  end subroutine countDiv

  subroutine countFunc( count )
    implicit none
    integer, intent(in) :: count

    if( itarget == 0 ) then
       write(*,*) "Error: count has not started yet."
       stop
    end if

    !$OMP CRITICAL
    vcount(itarget)%count_func &
         = vcount(itarget)%count_func + count
    !$OMP END CRITICAL

    return
  end subroutine countFunc

end module count_module

!!$program test
!!$  use count_module
!!$  implicit none
!!$
!!$  integer :: i, j, n
!!$  real(8) :: sqrt, exp, sin, erfc
!!$  real(8) :: a
!!$
!!$  n=1024
!!$
!!$  call initCount
!!$
!!$  call startCount("test1")
!!$  do i=1, n
!!$     call startCount("test2")
!!$     do j=1, n
!!$        a = sqrt(dble(i))*exp(-dble(j))+sin(dble(i+j))*erfc(dble(i-j))
!!$     end do
!!$     call countLoad(n*4)
!!$     call countStore(n)
!!$     call countAdd(n*1)
!!$     call countMul(n*2)
!!$     call countFunc(n*4)
!!$     call stopCount
!!$  end do
!!$  call countLoad(n*n*4)
!!$  call countStore(n*n)
!!$  call countAdd(n*n*1)
!!$  call countMul(n*n*2)
!!$  call countFunc(n*n*4)
!!$  call stopCount
!!$
!!$  call showCount
!!$
!!$  call exit(0)
!!$end program test


