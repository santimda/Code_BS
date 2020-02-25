
Program conv
  ! This program makes the convolution of the emission in the plane of the sky 
  ! with the instrument response 
  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(dp), parameter :: pi = acos(-1.d0)
  real(dp), parameter :: f_Ny = 0.6d0
  integer, parameter :: m_x_grilla = 50
  integer, parameter :: m_y_grilla = 50
  integer, parameter :: m_tot = m_x_grilla*m_y_grilla ! Total number of data
  real(dp), parameter :: res_x = 1.d0, res_y = 1.d0 ! Resolution 
  real(dp), parameter :: coc = 2.d0         
  real(dp), parameter :: theta_x = res_x*coc, theta_y = res_y*coc ! Resolution at other freq.
  real(dp) :: x_map_aux(m_x_grilla),y_map_aux(m_x_grilla)
  real(dp) :: x_map(m_tot),y_map(m_tot)
  real(dp) :: F_1Gabs(m_tot),D_x_grilla,D_y_grilla
  real(dp) :: ancho_x, ancho_y, D_x, D_y, sigma_x, sigma_y,ci
  real(dp) :: x_c0,y_c0,d_pesada2,suma,D_f,atenua,sum,mJy_beam,norm
  real(dp), dimension(:), allocatable :: x_c,y_c
  real(dp), dimension(:,:), allocatable :: F_map1Gabs
  integer :: m_x,m_y
  integer :: i, j, l
  real(dp) :: A_beam_elipse, A_beam_celda, A_grilla
  
  A_beam_elipse = f_Ny**2.*theta_x*theta_y
  A_beam_celda = (pi/4.d0)*theta_x*theta_y

  ancho_x = 20.d0  
  ancho_y = 20.d0 
  D_x_grilla = ancho_x/float(m_x_grilla)
  D_y_grilla = ancho_x/float(m_x_grilla)

  A_grilla = D_x_grilla*D_y_grilla 
  
  write(*,*) A_beam_elipse,A_beam_celda,A_grilla

  ! Región cuadrada de -10 a 10
  do i = 1,m_x_grilla
     x_map_aux(i) = -10.d0+(i-0.5)*D_x_grilla
!     write(*,*) x_map(i)
  enddo
  do j = 1,m_y_grilla
     y_map_aux(j) = -10.d0+(j-0.5)*D_y_grilla
!     write(*,*) y_map(j)
  enddo

  !2 emisores puntuales en el (-D/4,0),(D/4,0), con F = 100 c/u
  open(13,file='SW_radiomap_1GHz.dat')
  F_1Gabs = 0.d0
  l = 0
  do i = 1,m_x_grilla
     do j = 1,m_y_grilla
        l = l + 1
        x_map(l) = x_map_aux( mod(i,m_x_grilla)+1 )
        y_map(l) = y_map_aux( mod(j,m_y_grilla)+1 )
!        write(*,*) i,j,l,mod(i,m_x_grilla)+1,mod(j,m_y_grilla)+1
        !if ( l == (m_tot-m_x_grilla/2)/2 .or. l==(m_tot+m_x_grilla/2)/2 ) then    
!        if ( y_map(l) >= 0.d0 .and. abs(y_map(l)) < D_y_grilla .and. &
!             abs( abs(x_map(l)) - 2.5d0) < D_x_grilla .and. x_map(l) >= 0.d0 ) then
        if ( abs(y_map(l)) < 2.*D_y_grilla .and. (abs(x_map(l)) - 2.5d0) < D_x_grilla ) then
           !        if ( mod(j,m_y_grilla)+1 == m_y_grilla/2 ) then
           F_1Gabs(l) = 100.d0 !/A_grilla ! per beam
           write(*,*) l,x_map(l),y_map(l),F_1Gabs(l)
        endif
!        write(*,*) l,x_map(l),y_map(l)
        write(13,*) x_map(l),y_map(l),F_1Gabs(l)
     enddo
  enddo
  close(13)
!  write(*,*) l,x_map(l),y_map(j)

  !///////////////////////////////////////////////////////
  !-------------------------------------------------------
  ! Initialize the data for the convolution region
  
  
  ! theta is the HPBW, so we transform to sigma: (LOS MAPAS MUESTRAN SIGMA O THETA?!!!!!)
  sigma_x = theta_x/sqrt(8.d0*log(2.d0))
  sigma_y = theta_y/sqrt(8.d0*log(2.d0))
  norm = sigma_x*sigma_y*2.d0*pi
  
!  D_x = f_Ny*2.d0*sigma_x
!  D_y = f_Ny*2.d0*sigma_y
  D_x = f_Ny*theta_x
  D_y = f_Ny*theta_y
  m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
  m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
  write(*,*) m_x,m_y
  x_c0 = -ancho_x/2.d0 !+ D_x/2.d0  ! Center of the first ellipse 
  y_c0 = -ancho_y/2.d0 !+ D_y/2.d0
  allocate( x_c(m_x), y_c(m_y), F_map1Gabs(m_x,m_y)  )
  
  !//////////////////////////////////////////////////////
  !------------------------------------------------------
  ! Data read
  
  sum = 0.d0
  open(14,file='SW_radiomap_1GHz.dat')
  do l=1,m_tot
     read(14,*) x_map(l),y_map(l),F_1Gabs(l)
     sum = sum + F_1Gabs(l) !*A_grilla
  enddo
  write(*,*) 'Pre-conv:'
  write(*,80) 'S(1.4GHz)=',sum,'mJy'
  close(14)
    
  do i=1,m_x
     x_c(i) = x_c0 + (i-0.5)*D_x  ! x-coord of the center of the i-th ellipse
     
     do j=1,m_y
        y_c(j) = y_c0 + (j-0.5)*D_y  ! y-coord of the center of the (x_i,y_j) ellipse
        
        suma=1.d-9  ! simulate noise (?)
        do l=1,m_tot
           d_pesada2 = ((x_map(l)-x_c(i))**2.d0)/(2.d0*sigma_x**2.d0) + &
                ((y_map(l)-y_c(j))**2.d0)/(2.d0*sigma_y**2.d0) ! added extra sigma
           if (d_pesada2 < 6.d0) then      
              ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
              atenua = exp( -d_pesada2 ) / norm
              D_F = F_1Gabs(l)*atenua
           else
              D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
           endif
           suma = suma + D_f
        enddo
        F_map1Gabs(i,j) = suma
     enddo
     
  enddo
  
  !///////////////////////////////////////////////////////////////
  !
  ! Data storage + consistency check:
  
  ! Transform mJy --> mJy/beam by multiplying by A_beam/A_cell=pi/(2f_Ny)**2
  !    mJy_beam = pi/(2.*f_Ny)**2.d0 
!  mJy_beam = pi/f_Ny**2.d0 
  mJy_beam = A_grilla/A_beam_elipse
  
  sum = 0.d0
  open (17,file='SW_conv_1Gabs.dat')
  do i=1,m_x
     do j=1,m_y
        write(17,60) x_c(i),y_c(j),F_map1Gabs(i,j)*A_beam_elipse!*mJy_beam
        sum = sum + F_map1Gabs(i,j)*A_beam_elipse
     enddo
     write(17,*) 
  enddo
  write(*,*) 'Post-conv:'
  write(*,80) 'S(1.4GHz)=',sum,'mJy'
  write(*,*) ''
  close(17)
  
60 format(F11.4,X,F11.4,X,E11.4)
80 format(A11,X,Es10.2,X,A3)
  
!  write(*,*) m_x_grilla*m_y_grilla*D_x_grilla*D_y_grilla
!  write(*,*) float(m_x)*float(m_y) * f_Ny**2.*theta_x*theta_y

  deallocate (x_c,y_c,F_map1Gabs)
  
END Program conv
