module CD
  
  use global
  
  implicit none
  private
  public :: CD_run
  
contains
  
  !============================================================================
  ! CD_run
  !============================================================================
  subroutine CD_run
    
    !иииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииии
    ! This module calulates the 2D curve that represents the position of the shock  
    ! according to Christie+2016. Only the "upper" half is calculated, as there is 
    ! azimuthal symmetry. The "rotated" positions are later calculated by CDrot.f90
    !
    !иииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииииии
 
    use global
    implicit none
      

    open(11,file='CD.dat')

    do l=1,m_theta

       z(l) = R0*float(l)/float(m_theta)
       dR(l) = R0/float(m_theta)
       write(11,*) l,z(l),dR(l)

    enddo   
   
    close(11)
    
  End subroutine CD_run
  

End module CD

