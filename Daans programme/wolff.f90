PROGRAM wolff
    
! I am using KB=1 and J=1 then Beta*J = 1/T
  
  use initialconditions
  use cluster

  implicit none
  integer, parameter :: Niter=300, Tindex=500, DimX=100,DimY=100  ! number of iterations, X dimenssion, Y dimenssion 
  real(8) :: T = 0.01, BJ                                    ! Temperature

  integer :: M(DimY,DimX)=0, Maginit                         ! spin
  logical :: cluster(DimY,DimX)                              ! cluster wolff algorith
  real(8) :: addprobability                                  ! 1-exp(-2J/kT) probability to be part of the cluster
  integer :: i,j,x,y                                         ! extra variables
  !variable to measure chi and its error
  real(8) :: chi=0.0, chisum=0.0, chisqdsum=0.0 
  integer :: nchi
  
  !initialize
  call print(DimX,DimY,M(:,:),2)
  Maginit = sum(M(:,:))                                      ! calculate the initial magnetization 
  call init(DimY,DimX,M)                                     ! initialize ramdomly the spins 

  
  do j=1,Tindex
     T = 0.1+0.01*j
     BJ = 1.0/T                                              ! BJ factor, KB=1, J=1, then BJ= 1/T
     addprobability = 1-exp(-2.0*BJ)
     cluster(:,:)= .FALSE.
     do i=1,Niter
        x=nint( (DimX-1)*rand(0))+1
        y=nint( (DimY-1)*rand(0))+1
        call growcluster(x,DimX,y,DimY,M(:,:),cluster(:,:),addprobability)
     end do
     chi = (real(sum(M(:,:))))/((Dimx-1)*(Dimy-1))
     chisum = chisum + chi
     chisqdsum = chisqdsum + chi*chi
     nchi = nchi + 1
     print *,BJ,T, chi
  end do
  
  call print(DimX,DimY,M(:,:),0)

contains

subroutine print(DimX,DimY,M,i)
  integer, intent(in) :: DimX, DimY, M(:,:),i
1 FORMAT(' ',3I5)
  if(i>0) then
     open (unit = 4, file = 'initial.pnm')
     write (4,'(A)') "P1"
     write (4,1) DimX,DimY
!     write (4,*) "2"
     write (4,1) M(:,:)+1
  else
     open (unit = 5, file = 'final.pnm')
     write (5,'(A)') "P1"
     write (5,1) DimX,DimY
!     write (5,*) "2"
     write (5,1) M(:,:)+1
  end if
end subroutine print

END PROGRAM wolff
