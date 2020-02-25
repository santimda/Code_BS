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
    
    use global
    implicit none
    !¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    !     Este programa se encarga de calcular una curva 2D en la que está la 
    !     discontinuidad de contacto de una región de colisión de vientos. Tan
    !     sólo calcula la mitad superior, la inferior (y sus versiones rotadas)
    !     luego lo hace CDrot.f
    !________________________________________________________________________
    real(dp) :: dx_prom,dx,dl_max_ini,dl_max,dl,L_total
    allocate(x(ml),y(ml),z(ml))

    ! Definimos todos los parámetros del sistema(1 = O2If, 2 = O3.5V) sacados 
    ! de Muijres 2012. Para la O3.5 puedo interpolar, la otra requiere 
    ! extrapolar y es menos confiable, por eso intento no usar esos datos.
    
    !////////////////////////////////////////////////////////////////////////
    !
    !     Esquema: Busco los puntos (x,y) de la CD. Para ello voy a armar un 
    !     vector x(n),y(n), tal que (x(1),y(1)) sean las coordenadas del apex
    !     y los posteriores valores sean la rama superior. En (0,0) está la 
    !     estrella 1 y en (D,0) está la estrella 2. r1 y r2 son las distancias
    !     del punto (x,y) a la estrella 1 y 2 respectivamente. Trabajaré el
    !     caso simplificado en que los vientos ya alcanzaron su velocidad 
    !     terminal y por tanto eta=cte. Desprecio movimiento orbital por lo que
    !     la rama inferior será el espejo de la inferior.
    
    !      y(1)=0.d0 --> Si uso esto tengo problemas numéricos
    
    y(1)=1.d-5
    x(1)=(D/UA)/(1.d0+sqrt(eta))
    ! Del global lee: x_max prop to (D/UA)
    dx_prom = 2.d0*(x_max-x(1))/float(ml/2)
    dl_max_ini = 9.d11/UA/cos(incli)
    ! Es la que utilizo para el primer paso, luego agrando
    ! un poco ya que no necesito tanta resolución lejos del apex
    L_total = 0.d0
    do i=2,ml/2
       !Agrando un poco el paso
       dl_max = dl_max_ini*(1.d0+30.d0*float(i-2)/float(ml))
       x(i)=x(i-1)+dx_prom
       !     La subrutina RK me me permite obtener y(i) pasándole como argumento
       !     x(i-1),y(i-1) y dx. La función dy/dx ya la tiene como externa (g), 
       !     pero es necesario actualizarle a mano D, eta en g.
       dl = dl_max + 1.d0 ! me aseguro que entre
       do while (dl > dl_max)
          dx = (x(i)-x(i-1))/2
          x(i)= x(i-1) + dx
          call RK(x(i-1),y(i-1),dx,y(i))
          dl = sqrt(dx**2+(y(i)-y(i-1))**2)
       enddo
       dx_prom = 2.d0*(x_max-x(i))/float((ml/2)-i) ! Re-adapto el paso
!       write(*,*) dl_max/dl_max_ini, dl/dl_max_ini
       L_total = L_total + dl
    enddo
    write(*,*) 'L_total [D] =', 2.d0*L_total/(D/UA)

    open(11,file='CD.dat')
    do l=1,ml/2
       write(11,*) x(l)*UA,y(l)*UA ! La coord z es cero, luego roto
    enddo
    do l=ml/2+1,ml
       write(11,*) x(l-ml/2)*UA,-y(l-ml/2)*UA ! La coord z es cero, luego roto
    enddo
    close(11)
    
    deallocate (x,y,z)
  End subroutine CD_run
  

End module CD



!/////////////////////////////////////////////////////////////////

Subroutine RK(xi,yi,dx,yf)
  
  implicit none
  
  integer, parameter :: dp = kind(1.d0)
  real(dp), intent(in) :: xi,yi,dx
  real(dp), intent(out) :: yf
  real(dp) :: y,k1,k2,k3,k4,h,g
  integer :: i
  integer, parameter :: n=4000
  dimension y(n)
  external g
  
  y(1)=yi
  h=dx/float(n)
  
  do i=1,n
     k1=h*g(xi+(i-1)*h,y(i))
     k2=h*g(xi+(i-1)*h+h/2.d0,y(i)+k1/2.d0)
     k3=h*g(xi+(i-1)*h+h/2.d0,y(i)+k2/2.d0)
     k4=h*g(xi+(i)*h,y(i)+k3)
     y(i+1)=y(i)+(k1+2.d0*k2+2.d0*k3+k4)/6.d0
  enddo
  
  yf=y(n)
  
  return
  
End subroutine RK

!/////////////////////////////////////////////////////////////////


Function g(a_i,a_f)
  
  use global
  
  real(dp) :: g
  real(dp) :: r1,r2,a_i,a_f
  r1=sqrt( a_i**2 + a_f**2 )
  r2=sqrt( (D/UA - a_i)**2 + a_f**2 )
  g=a_f/(a_i-(D/UA*r1**2*sqrt(eta)/(r1**2*sqrt(eta)+r2**2)))! dy/dx
  return
End function g

!/////////////////////////////////////////////////////////////////
