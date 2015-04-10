module susceptibility

    implicit none
    private

    public chimeasure

contains
  
  subroutine chimeasure(DimY,DimX,M)
    integer, intent(in) :: DimY,DimX
    integer, intent(inout) :: M(:,:)
    real(8) :: chi, chisum, chisqdsum
    integer :: nchi
    
    

  end subroutine chimeasure
  
end module
