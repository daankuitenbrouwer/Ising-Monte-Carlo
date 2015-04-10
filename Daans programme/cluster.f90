module cluster

  implicit none
  private
  
  public growcluster, tryadd
  
contains
  
  subroutine growcluster(x,DimX,y,DimY,M,cluster,addprobability)
    integer, intent(in) :: DimX,DimY
    integer, intent(inout) :: x,y,M(:,:)
    logical, intent(inout) :: cluster(:,:)
    real(8), intent(in) :: addprobability
    integer :: xprev, xnext, yprev, ynext,currentspin
       
    cluster(y,x) = .TRUE.
    currentspin = M(y,x)
    M(y,x) = -M(y,x)
     
    xprev = x-1
    xnext = x+1
    yprev = y-1
    ynext = y+1
    
    if( xprev == 0 ) xprev = DimX
    if( xnext == DimX+1) xnext = 1
    if( yprev == 0 ) yprev = DimY
    if( ynext == DimY+1) ynext = 1

    if(cluster(y,xprev).eqv. .FALSE.)then
       call tryadd(xprev,DimX,y,DimY,M,cluster,addprobability,currentspin)
    end if

    if(cluster(y,xnext).eqv. .FALSE.)then
       call tryadd(xnext,DimX,y,DimY,M,cluster,addprobability,currentspin)
    end if

    if(cluster(yprev,x).eqv. .FALSE.)then
       call tryadd(x,DimX,yprev,DimY,M,cluster,addprobability,currentspin)
    end if

    if(cluster(ynext,x).eqv. .FALSE.)then
       call tryadd(x,DimX,ynext,DimY,M,cluster,addprobability,currentspin)
    end if

  end subroutine growcluster
  

  subroutine tryadd(x,DimX,y,DimY,M,cluster,addprobability,currentspin)
    integer, intent(in) :: DimX,DimY,currentspin
    integer, intent(inout) :: x,y,M(:,:)
    real(8), intent(in) :: addprobability
    logical, intent(inout) :: cluster(:,:)
    
    if (M(y,x)==currentspin) then
       if(rand(0)<addprobability) then
          call growcluster(x,DimX,y,DimY,M,cluster,addprobability)
       end if
    end if
    
  end subroutine tryadd
  

End module cluster
