program FlatSky 
  implicit none
  integer, parameter :: dl= KIND(1.d0)
  real(dl), parameter :: pi = 3.14159265359

  character(80) :: Folder1, Folder2, Folder3
  character(80) :: Clfile, Cllfile

  real(dl), pointer :: Cl(:,:), Cll(:,:)
  real(dl), pointer :: pClpp(:,:)
  integer :: l1, l2a, l3a, l2b, l3b
  integer :: l3blmin, l3blmax,l3almax,l3almin, l1a,l1b
  integer :: lmax, lmin
  integer :: i,j,k,m
  real(dl) :: phi23a, phi23b, phi31a, phi31b, phi21a, phi21b, phi2a3b, phi2b3a, phi2a2b, phi3a3b
  real(dl) :: fnl, signsq
  real(dl) :: absl2l3a, absl2l3b, absl3al3b, absl2al2b, absl2al3b, absl2bl3a
  real(dl) :: l2dotl3a, l2dotl3b, l2adotl3b,l2bdotl3a,l3adotl3b,l2adotl2b,l1dotl2a,l1adotl2b
  real(dl) :: DB(3)
  real(dl) :: SumGauss, DSNGauss,  SumNGauss, DSNonGauss
  real(dl) :: CMB2COBEnorm = 7428350250000.d0

  integer :: ellar(512), dellar(512), elltoi(5000)
  integer :: intmax, imin
  integer :: stp



  !various binning schemes
!!$  do i  = 1, 256
!!$     if (i .le. 19) then
!!$        ellar(i) = i+1
!!$        dellar(i) = 1.d0
!!$     elseif (i .le. 75 .and. i .ge. 20) then 
!!$        ellar(i) = ellar(i-1) + 4
!!$        dellar(i) = 4.d0
!!$     elseif (i .le. 257 .and. i .ge. 76) then
!!$        ellar(i)  = ellar(i-1) + 20
!!$        dellar(i) = 20.d0
!!$     endif
!!$     !write(*,*) 'ell:', ellar(i)
!!$  enddo
  do i  = 1, 256
     if (i .le. 47) then
        ellar(i) = i+1
        dellar(i) = 1.d0
        elltoi(ellar(i)) = i
     elseif (i .le. 85 .and. i .ge. 48) then 
        ellar(i) = ellar(i-1) + 4
        dellar(i) = 4.d0
        elltoi(ellar(i)) = i

     elseif (i .le. 110 .and. i .ge. 86) then
        ellar(i)  = ellar(i-1) + 12
        dellar(i) = 12.d0
        elltoi(ellar(i)) = i

     elseif (i .le. 175 .and. i .ge. 111) then
        ellar(i)  = ellar(i-1) + 24
        dellar(i) = 24.d0
        elltoi(ellar(i)) = i

     else
        ellar(i)  = ellar(i-1) + 50
        dellar(i) = 50.d0
        elltoi(ellar(i)) = i

     endif
     !write(*,*) 'ell:', ellar(i)
  enddo
!!$  do i  = 1, 512
!!$     if (i .le. 399) then
!!$        ellar(i) = i + 1
!!$        dellar(i) = 1.d0
!!$     else
!!$        ellar(i)  = ellar(i-1) + 2
!!$        dellar(i) = 2.d0
!!$     endif
!!$     !write(*,*) 'ell:', ellar(i)
!!$  enddo
  !stop
  !∆`=  1  for`≤50,  ∆`=  4  for50< `≤200,  ∆`=  12  for  200< `≤500,  ∆`=  24for 500< `≤2000,  and finally ∆`= 40 for` >2000   

  lmax = 5000

  lmin = 2

  allocate(Cl(4,2:lmax))
  allocate(Cll(4,2:lmax))
  allocate(pClpp(3,2:lmax))

  Folder1 = './SOspectra/'
  Clfile = trim(Folder1)//trim('SOspectra_lenspotentialCls.dat')
  Cllfile = trim(Folder1)//trim('SOspectra_lensedCls.dat')
  !from Alex:
  Cllfile = './SO_forecasts/CAMB/cosmo2017_10K_acc3_lensedCls.dat'
  open(unit=17,file = Clfile, status='old')
  open(unit=18,file = Cllfile, status='old')

  do j = 1, lmax
     !#    L    TT             EE             BB             TE 
     !#    L    TT             EE             BB             TE             PP             TP             EP
     if (j .eq. 1) then
        read(17,*)
        read(18,*)
        cycle
     endif

     read(17,*) l1, Cl(1:4,j),pClpp(1:3,j)
     read(18,*) l1, Cll(1:4,j)
     Cll(1:4,j) = 2.*pi*Cll(1:4,j)/l1/(l1+1.)/CMB2COBEnorm
     Cl(1:4,j) = 2.*pi*Cl(1:4,j)/l1/(l1+1.)/CMB2COBEnorm
     pClpp(1,j) = 2.*pi*pClpp(1,j)/(l1*(l1+1.))**2
     pClpp(2:3,j) = 2.*pi*pClpp(2:3,j)/(l1*(l1+1.))**(3.d0/2.d0)/CMB2COBEnorm**(1./2)

     !write(*,*) l1,Cll(1:4,j) 
  enddo
  close(17)
  close(18)
  SumGauss =0.d0
  SumNGauss = 0.d0
  DB = 0.d0
  !structure
  !l2/l3/l3/l2'/l3'
  lmax = 400
  lmin = 2

  !intmax = 60 : lmax = 100
  intmax = 60

  imin = 9
  lmax = ellar(intmax)
  lmin = ellar(imin)

  !to make it faster:
  stp  = 1

  write(*,*) 'lmax:',lmax
  write(*,*) 'lmin:',lmin

  j = 1
  

  !do i = imin, intmax !l2 loop
  do i = imin, intmax   
     l1a = ellar(i)

     !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(dynamic) &
     !$OMP PRIVATE(l1b,l2a,l3a,l2b,l3b,l3almax,l3almin,l3blmin,l3blmax,j,k), &
     !$OMP PRIVATE(fnl,phi23a,phi23b,signsq,phi21a,phi31a,phi21b,phi31b,phi2a3b,phi2b3a), &
     !$OMP PRIVATE(phi2a2b,phi3a3b,l2dotl3a,l2dotl3b,l2adotl3b,l2bdotl3a,l3adotl3b,l2adotl2b), &
     !$OMP PRIVATE(absl2al3b,absl2al2b,absl3al3b,absl2bl3a,DB,DSNonGauss,DSNGauss,l1dotl2a,l1adotl2b), &
     !$OMP REDUCTION(+:SumGauss,SumNGauss)
     do j = imin, intmax !l3 loop
     !do l3a = lmin, lmax   
        l2a = ellar(j)
        !set minimum and max values of third leg (triangle constraint)
        l3almax = min(l2a+l1a,lmax)
        l3almin = max(abs(l2a-l1a),lmin)

        do l3a = l3almin, l3almax, stp

           !angle between l2 and l3
           phi23a = angle(l2a,l3a,l1a)
           !inner products:
           !l2dotl3a = l2a*l3a*cos(phi23a) !l2.l3
           l1dotl2a = l1a*l2a*cos(phi21a)

           !other angles; needed for 11 permutations 
           phi21a = angle(l1a,l2a,l3a)
           phi31a = pi - phi23a - phi21a

           !get the SW bispectrum for the triplet (l1, l2 ,l3)
           fnl = floc(l1a, l2a,l3a)

           !Eauate the l1 and l1'
           l1b = l1a
           !l2b can max be size of l1 (or l1'). Can be small, given l3
           !il2max =
           do k = imin, intmax
              l2b = ellar(k)
           !do k = imin, intmax !l2' loop
           !   l2b = ellar(k)

              !minimum value is contraint to
              l3blmin = max(abs(l1b-l2b),lmin)
              !maximum:
              l3blmax = min(l1b+l2b,lmax)

              do l3b = l3blmin, l3blmax, stp !final loop: l3'

                 !prime triangle angle between l2' and l3'
                 phi23b = angle(l2b,l3b,l1b)
                 !compute fnl^2
                 signsq = fnl*floc(l1b,l2b,l3b)

                 !other angles; needed for 11 permutations 
                 phi21b = angle(l1b,l2b,l3b)
                 l1adotl2b = l1a*l2b*cos(phi21b)
                 phi31b = pi - phi23b - phi21b

                 phi2a3b = pi - phi21a + phi31b !correct direction
                 phi2b3a = pi - phi21b + phi31a

                 phi2a2b = pi- phi21a + phi21b
                 phi3a3b = pi- phi31a + phi31b

                 !other inner products:
                 l2dotl3b = l2b*l3b*cos(phi23b) !l2'.l3'

                 l2adotl3b = l2a*l3b*cos(phi2a3b) !l2.l3'
                 l2bdotl3a = l2b*l3a*cos(phi2b3a) !l2'.l3

                 l3adotl3b = l3a*l3b*cos(phi3a3b) !l3.l3'
                 l2adotl2b = l2a*l2b*cos(phi2a2b) !l2.l2'

                 !other lengths:
                 absl2al2b = nint(length(l2a,l2b,phi2a2b)) !|l2+l2'|
                 absl3al3b = nint(length(l3a,l3b,phi3a3b)) !|l3+l3'|                 

                 absl2al3b = nint(length(l2a,l3b,phi2a3b)) !|l2+l3'|
                 absl2bl3a = nint(length(l2b,l3a,phi2b3a)) !|l2'+l3|
!!$                 
!!$                 
!!$
                 !if(absl2bl3a .lt. lmin) absl2bl3a = lmin
                 !if(absl2al2b .lt. lmin) absl2al2b = lmin

                 !if(absl2bl3a .gt. lmax) absl2bl3a = lmax
                 !if(absl2al2b .gt. lmax) absl2al2b = lmax
!!$                 
                 !3 unique terms 
                 DB(1) =  Cl(1,l1a)*Cl(1,l2a)*Cl(1,l2b)*pClpp(1,l1a)* &
                                !(l2a**2+l2dotl3a)*(l2b**2+l2dotl3b)
                      (l1adotl2b*l1dotl2a)
                 !write(*,*) l2a, DB(1)
!!$
!!$                 
!!$                 if(absl2bl3a .lt. lmin .or. absl2bl3a .gt. lmax) then
!!$                    DB(2) = 0.d0
!!$                 else                    
!!$                    DB(2) = ( Cl(1,l1a)*Cl(1,l2a)*Cl(1,l2b)*pClpp(1,absl2bl3a)* &
!!$                         (l2a**2+l2adotl3b)*(l2b**2+l2bdotl3a))
!!$                 endif
!!$                 if(absl2al2b .lt. lmin .or. absl2al2b .gt. lmax) then
!!$                    DB(3) = 0.d0
!!$                 else                    
!!$                    DB(3) = (Cl(1,l1a)*Cl(1,l2b)*Cl(1,l3b)*pClpp(1,absl2al2b)* &
!!$                         (l3b**2+l3adotl3b)*(l2b**2+l2adotl2b))       
!!$                 endif
!!$                          
!!$
                 DSNonGauss = signsq*Sum(DB(1:1))/Cl(1,l1a)**2/Cl(1,l2a)/Cl(1,l2b)/Cl(1,l3a)/Cl(1,l3b)/(2.*pi)**4 *dellar(i)*dellar(j)*(stp)**3
!!$
!!$                 !if (((l2a.eq.l2b) .and. (l3a .eq.l3b)) .or. ((l2a.eq.l3b) .and. (l3a .eq.l2b))) then
                 if (((l2a.eq.l2b) .and. (l3a .eq.l3b))) then
                    DSNGauss = signsq/Cl(1,l1a)/Cl(1,l2a)/Cl(1,l3a)/6./(2.*pi)**2 *dellar(i)*dellar(j)*stp
                    !<N^2> + delta <N^2>
                 else
                    DSNGauss = 0.d0
                    !<N^2> + delta <N^2>

                 endif
                 SumGauss = SumGauss + DSNGauss
                 SumNGauss =  SumNGauss + DSNGauss + DSNonGauss

              enddo !l3b

           enddo !l2b
        enddo !l1a
     enddo !l3a
     !$OMP END PARAllEl DO
     write(*,*) 'l2 = ', l1a
  enddo !l2a

  write(*,'(I4,4E17.8)') lmax, SumGauss, SumNGauss, sqrt(SumGauss/SumNGauss)
  write(*,'(A12,X,I4,X,A19,X,F11.3)') 'For lmax = ',  lmax, 'the error on fnl = ', sqrt(1./SumGauss)

  close(12)
  deallocate(Cl, Cll, pClpp)
contains

  real(dl) function floc(l1,l2,l3)
    !SW approximation 
    integer :: l1, l2, l3
    real(dl) :: amp
    real(dl) :: As = 2.1056d-9
    !from https://arxiv.org/pdf/0812.3413.pdf Eq. 19 and 20
    amp = (2.d0/27./pi**2)*As
    floc = 1.d0/(l1+1.d0)/l1/l2/(l2+1.d0) + 1.d0/(l3+1.d0)/l3/l2/(l2+1.d0) + &
         1.d0/(l1+1.d0)/l1/l3/(l3+1.d0)
    floc = floc*amp*2.E-7 !2.E-7  is introduced to get roughly same amplitude at l_max = 500 to full fnl_local
  end function floc

  real(dl) function length(x,y,phi)
    integer, intent(in) :: x, y
    real(dl), intent(in) :: phi
    length = sqrt(x**2+y**2-2.*x*y*cos(phi))
  end function length

  real(dl) function angle(x,y,z)
    integer, intent(in) :: x,y,z
    angle = acos((x**2+y**2-z**2)/(2.*x*y))
    !if (angle .ne. angle) angle  = 0.d0
  end function angle


end program FlatSky
