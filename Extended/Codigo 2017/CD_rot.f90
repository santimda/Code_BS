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
  
  subroutine CD_rot_run  
    use global
    implicit none
    real(dp) :: xrot,yrot,zrot,phi
    allocate ( x(ml),y(ml),z(ml) )
   
    open(11,file='phi.dat')
    read(11,*) phi            ! Está en rad
    close(11)

    open(12,file='CD.dat')
    open(13,file='CDrot.dat')
    do l=1,ml/2
       read(12,*) x(l),y(l)
       z(l)=0.d0
       call rotacion(x(l),y(l),z(l),phi,xrot,yrot,zrot)
       write(13,*) xrot,yrot,zrot
    enddo
    do l=ml/2+1,ml
       call rotacion(x(l-ml/2),-y(l-ml/2),0.d0,phi,xrot,yrot,zrot)
       write(13,*) xrot,yrot,zrot
    enddo
        close(12)
    close(13)

    deallocate ( x,y,z )
    
  END subroutine CD_rot_run
  
END module CD_rot
  
!____________________________________
!____________________________________
  
  
  Subroutine rotacion(Vx,Vy,Vz,phi,VRx,VRy,VRz)
    ! Rota el vector V un ángulo phi en torno al eje x
    real*8 vx,vy,vz,phi,vrx,vry,vrz
      
    VRx = vx
    VRy = cos(phi)*Vy ! Vz=0
    VRz = sin(phi)*Vy ! Vz=0
    
    RETURN
  END Subroutine rotacion
    
