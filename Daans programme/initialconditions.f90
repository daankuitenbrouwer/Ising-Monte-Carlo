module initialconditions

    implicit none
    private

    public init

contains
  
  subroutine init(DimY,DimX,M)
    integer, intent(in) :: DimY,DimX
    integer, intent(inout) :: M(:,:)
    integer :: ii,jj
    integer :: random

!    M(:,:)=-1
    do jj=1,DimY
       do ii=1,DimX
          random = -1!2*nint(rand(0))-1
          M(jj,ii)=random
       end do
    end do

  end subroutine init
  
end module
