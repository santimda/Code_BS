module CD_rot
  !Este módulo rota los vectores posición de la CD en torno al eje que 
  !une las 2 estrellas (eje x) calcula la emisión IC ya corregida por absorción
  use global
  
  implicit none
  private
  public :: CD_rot_run

contains

  !============================================================================
  ! CD_rot_run
  !============================================================================
  
  subroutine CD_rot_run(phi)  
    use global
    implicit none
    real(dp), intent(in) :: phi
    real(dp) :: xrot,yrot,zrot
    
    interface
       subroutine rotacion(x1,x2,x3,ang,x1_rot,x2_rot,x3_rot)
         use iso_fortran_env, only: dp => real64
         implicit none
         real(dp), intent(in) :: x1,x2,x3,ang
         real(dp), intent(out) :: x1_rot,x2_rot,x3_rot
       end subroutine rotacion
    end interface
    
    allocate( x(m_theta),y(m_theta),z(m_theta) )
              
     write(*,*) 'Rotate the coords'

    open(12,file='CD_xy.dat')
    open(13,file='CDrot.dat')
    do l = 1,m_theta
       read(12,*) x(l),y(l)
       z(l)=0.d0
       call rotacion(x(l),y(l),z(l),phi,xrot,yrot,zrot)
       write(13,*) xrot,yrot,zrot
    enddo
    do l = m_theta+1, 2*m_theta
       call rotacion(x(l-m_theta),-y(l-m_theta),0.d0,phi,xrot,yrot,zrot)
       write(13,*) xrot,yrot,zrot
    enddo
    close(12)
    close(13)
    
    deallocate( x,y,z )
    
  END subroutine CD_rot_run
  
END module CD_rot

!____________________________________
!____________________________________


Subroutine rotacion(Vx,Vy,Vz,ang,VRx,VRy,VRz)
  use global
  implicit none 
  ! Rotates a vector V in an angle "ang" around the x-axis
  real(dp), intent(in) :: Vx,Vy,Vz,ang
  real(dp), intent(out) :: VRx,VRy,VRz
      
  VRx = Vx
  VRy = cos(ang)*Vy ! Vz=0
  VRz = sin(ang)*Vy ! Vz=0
  
  RETURN
END Subroutine rotacion

