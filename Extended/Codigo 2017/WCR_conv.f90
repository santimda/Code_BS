module WCR_conv
  ! Este programa convoluciona la emisión de cada punto con la respuesta del 
  ! instrumento
  use global
  
  implicit none
  private
  public :: WCR_conv2G_run,WCR_conv8G_run
  
contains
  
  !============================================================================
  ! WCR_conv2G_run
  !============================================================================
  
  subroutine WCR_conv2G_run  
    use global
    implicit none
    integer, parameter :: ml_tot = 2*(mphi-1)*ml   ! Cantidad total de datos
    real(dp), parameter :: f_Ny = 0.53d0      ! Frecuencia de muestreo
    real(dp), parameter :: theta_x = 15.d0, theta_y = 11.d0
!    real(dp), parameter :: theta_x = 31.d0, theta_y = 23.d0
    real(dp) :: x_map(ml_tot),y_map(ml_tot)
    real(dp) :: F_2Gabs(ml_tot),F_2Gun(ml_tot)
    real(dp) :: ancho_x, ancho_y, D_x, D_y, sigma_x, sigma_y
    real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
    real(dp), dimension(:), allocatable :: x_c,y_c
    real(dp), dimension(:,:), allocatable :: F_map2Gabs
    integer :: m_x,m_y


    !///////////////////////////////////////////////////////
    !-------------------------------------------------------
    ! Inicializo los datos de la región donde convolucionaré

    ancho_x = 4.d0*D_proj ! [mas]
    ancho_y = 4.d0*D_proj ! [mas]

    ! theta es el HPBW, y yo necesito sigma
    sigma_x = theta_x/sqrt(8.d0*log(2.d0))
    sigma_y = theta_y/sqrt(8.d0*log(2.d0))

    D_x = f_Ny*2.d0*sigma_x
    D_y = f_Ny*2.d0*sigma_y
    m_x = int(ancho_x/D_x)+1
    m_y = int(ancho_y/D_y)+1
    x_c0 = -ancho_x/2.d0  ! Donde centro la primera elipse
    y_c0 = -ancho_y/4.d0
    allocate( x_c(m_x+1), y_c(m_y+1), F_map2Gabs(m_x+1,m_y+1)  )


    !//////////////////////////////////////////////////////
    !------------------------------------------------------
    ! Ahora leo los datos

    sum = 0.d0
    open(13,file='WCR_radiomap_2GHz.dat')
    do l=1,ml_tot
       read(13,*) y_map(l),x_map(l),F_2Gun(l),F_2Gabs(l)
       sum = sum + F_2Gabs(l)
    enddo
    write(*,*) 'Pre-conv: S(2.3GHz)=',sum,'mJy'
    close(13)
    
 
    do i=1,m_x+1
       x_c(i) = x_c0 + (i-1)*D_x
       
       do j=1,m_y+1
          y_c(j) = y_c0 + (j-1)*D_y
          
          suma=1.d-9  ! simulo el ruido 
          do l=1,ml_tot
             d_pesada2 = ((x_map(l)-x_c(i))**2.d0)/(2.d0*sigma_x) + &
                  ((y_map(l)-y_c(j))**2.d0)/(2.d0*sigma_y)
             if (d_pesada2 < 8.d0) then              
                atenua = exp( -d_pesada2 )!/(sigma*sqrt(2*pi)) ?!
                D_F = F_2Gabs(l)*atenua
             else
                D_F = 0.d0
             endif
             suma = suma + D_f
          enddo
          F_map2Gabs(i,j) = suma
       enddo

    enddo

    !///////////////////////////////////////////////////////////////
    !
    ! Y ahora a escribir los datos (+ chequeo de consistencia):

    !Para pasar a mJy/beam debo multiplicar por A_beam/A_celda=pi/(2f_Ny)**2
    mJy_beam = pi/(2.d0*f_Ny)**2.d0 
    sum = 0.d0
    open (17,file='WCR_conv_2Gabs.dat')
    do i=1,m_x+1
       do j=1,m_y+1
          write(17,60) x_c(i),y_c(j),F_map2Gabs(i,j)*mJy_beam
          sum = sum + F_map2Gabs(i,j)
       enddo
       write(17,*) 
    enddo
    write(*,*) 'S(2.3GHz)=',sum,'mJy'
    close(17)
    
60  format(F11.4,X,F11.4,X,E11.4)
    
    deallocate (x_c,y_c,F_map2Gabs)

  END subroutine WCR_conv2G_run

  !============================================================================
  ! WCR_conv2G_run (END)
  !============================================================================
  
  
  !============================================================================
  ! WCR_conv8G_run
  !============================================================================
  
  subroutine WCR_conv8G_run  
    use global
    implicit none
    integer, parameter :: ml_tot = 2*(mphi-1)*ml ! Cantidad total de datos
    real(dp), parameter :: f_Ny = 0.99d0          ! Frecuencia de muestreo
    real(dp), parameter :: coc = 2.3d0/8.6d0     ! Resolución extra 
    real(dp), parameter :: theta_x = 15.d0*coc, theta_y = 11.d0*coc
    real(dp) :: x_map(ml_tot),y_map(ml_tot)
    real(dp) :: F_8Gabs(ml_tot),F_8Gun(ml_tot)
    real(dp) :: ancho_x, ancho_y, D_x, D_y, sigma_x, sigma_y
    real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam
    real(dp), dimension(:), allocatable :: x_c,y_c
    real(dp), dimension(:,:), allocatable :: F_map8Gabs
    integer :: m_x,m_y


    !///////////////////////////////////////////////////////
    !-------------------------------------------------------
    ! Inicializo los datos de la región donde convolucionaré

    ancho_x = 4.d0*D_proj ! [mas]
    ancho_y = 4.d0*D_proj ! [mas]

    ! theta es el HPBW, y yo necesito sigma
    sigma_x = theta_x/sqrt(8.d0*log(2.d0))
    sigma_y = theta_y/sqrt(8.d0*log(2.d0))

    D_x = f_Ny*2.d0*sigma_x
    D_y = f_Ny*2.d0*sigma_y
    m_x = int(ancho_x/D_x)+1
    m_y = int(ancho_y/D_y)+1
    x_c0 = -ancho_x/2.d0  ! Donde centro la primera elipse
    y_c0 = -ancho_y/4.d0
    allocate( x_c(m_x+1), y_c(m_y+1), F_map8Gabs(m_x+1,m_y+1)  )


    !//////////////////////////////////////////////////////
    !------------------------------------------------------
    ! Ahora leo los datos

    sum = 0.d0
    open(13,file='WCR_radiomap_8GHz.dat')
    do l=1,ml_tot
       read(13,*) y_map(l),x_map(l),F_8Gun(l),F_8Gabs(l)
       sum = sum + F_8Gabs(l)
    enddo
    write(*,*) 'Pre-conv: S(8.6GHz)=',sum,'mJy'
    close(13)
    
 
    do i=1,m_x+1
       x_c(i) = x_c0 + (i-1)*D_x
       
       do j=1,m_y+1
          y_c(j) = y_c0 + (j-1)*D_y
          
          suma=1.d-9  ! simulo el ruido 
          do l=1,ml_tot
             d_pesada2 = ((x_map(l)-x_c(i))**2.d0)/(2.d0*sigma_x) + &
                  ((y_map(l)-y_c(j))**2.d0)/(2.d0*sigma_y)
             if (d_pesada2 < 8.d0) then              
                atenua = exp( -d_pesada2 )
                D_F = F_8Gabs(l)*atenua
             else
                D_F = 0.d0
             endif
             suma = suma + D_f
          enddo
          F_map8Gabs(i,j) = suma
       enddo

    enddo

    !///////////////////////////////////////////////////////////////
    !
    ! Y ahora a escribir los datos + chequeo de consistencia:
    
    !Para pasar a mJy/beam debo multiplicar por A_beam/A_celda = pi/f_Ny**2
    mJy_beam = pi/f_Ny**2.d0 

    sum = 0.d0
    open (18,file='WCR_conv_8Gabs.dat')
    do i=1,m_x+1
       do j=1,m_y+1
          write(18,61) x_c(i),y_c(j),F_map8Gabs(i,j)*mJy_beam
         sum = sum + F_map8Gabs(i,j)
       enddo
       write(18,*) 
    enddo
    write(*,*) 'S(8.6GHz)=',sum,'mJy'
    close(18)

61  format(F11.4,X,F11.4,X,E11.4)
    deallocate (x_c,y_c,F_map8Gabs)

  END subroutine WCR_conv8G_run

  !============================================================================
  ! WCR_conv8G_run (END)
  !============================================================================
  
END module WCR_conv
