program FlatSky 
  implicit none
  integer, parameter :: dl= KIND(1.d0)
  type :: bfs
     real(dl), allocatable :: rarray(:), deltar(:),Cll(:,:),pCll(:,:)
     real(dl), allocatable :: ar(:,:,:), br(:,:,:), gr(:,:,:), dr(:,:,:)
     integer :: flagDoWigner
  end type bfs


  real(dl), parameter :: pi = 3.14159265359



  character(120) :: Folder1, Folder2, Folder3
  character(120) :: Clfile, Cllfile,suffix

  real(dl), pointer :: Cl(:,:), Cll(:,:), invCll(:,:,:),Cllm(:,:,:)

  real(dl), pointer :: pClpp(:,:)
  integer :: l1, l2,l3, l2a, l3a, l2b, l3b,max_l,min_l
  integer :: l3blmin, l3blmax,l3almax,l3almin, l1a,l1b
  integer :: lmax, lmin,lmaxTmp 
  integer :: i,j,k,m,n,m2,p2,q2,q1,p1,m1
  real(dl) :: phi23a, phi23b, phi31a, phi31b, phi21a, phi21b, phi2a3b, phi2b3a, phi2a2b, phi3a3b
  real(dl) :: phi21amax,phi21bmax,phi21amin,phi21bmin,dPhia,dPhib
  real(dl) :: fnl,fnlb, signsq
  real(dl) :: absl2l3a, absl2l3b, absl3al3b, absl2al2b, absl2al3b, absl2bl3a
  real(dl) :: Rl2adotl3b,Rl2bdotl3a,Rl3adotl3b ,Rl2adotl2b,Rabsl2al2b ,Rabsl2al3b
  real(dl) :: l2dotl3a, l2dotl3b, l2adotl3b,l2bdotl3a,l3adotl3b,l2adotl2b,l1dotl2a,l1adotl2b
  real(dl) :: DB(3),DB_tmp
  real(dl) :: SumGauss, DSNGauss,  SumNGauss, DSNonGauss, SumNGaussAll, measureG, measureNG
  real(dl) :: CMB2COBEnorm = 7428350250000.d0


  real(dl)  :: atj(0:20000),atj2(0:20000)
  real(dl), pointer :: a3j(:,:), a3joC(:,:)
  real(dl), pointer :: bispectrum(:,:,:,:,:),bispectrum_ISWlens(:,:,:,:,:), bis(:,:),bis_ISWlens(:,:)
  real(dl) ::  tempfac,tempfacFcM(2,2),tmpPrefac
  real(dl) ::  DetFishCV,TotSumCVISWLens,TotSumCV,TotSumCVLensCross,alpha,beta
  real(dl) ::  Det,DetISWLens,fnlISW,fnlISWb,TempCovCV,DetLensCross,detCovCV

  type(bfs) :: P
  type(bfs) :: P_ISWlens
  integer :: ellar(2048), dellar(2048), elltoi(5000)
  integer :: intmax, imin,nPhib,nPhia,nPhi
  integer :: stp

  integer ::  shape,nfields,minfields
  character(120) :: alphabetafile, alphabetaPolfile,tensDir

  logical :: want_ISW_correction = .false.
  logical :: want_all_terms = .true.

  shape = 1
  minfields = 1
  nfields = 1
  want_ISW_correction = .false.

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
  ! do i  = 1, 512
  !    if (i .le. 47) then
  !       ellar(i) = i+1
  !       dellar(i) = 1.d0
  !       elltoi(ellar(i)) = i
  !    elseif (i .le. 85 .and. i .ge. 48) then 
  !       ellar(i) = ellar(i-1) + 4
  !       dellar(i) = 4.d0
  !       elltoi(ellar(i)) = i

  !    elseif (i .le. 110 .and. i .ge. 86) then
  !       ellar(i)  = ellar(i-1) + 12
  !       dellar(i) = 12.d0
  !       elltoi(ellar(i)) = i

  !    elseif (i .le. 175 .and. i .ge. 111) then
  !       ellar(i)  = ellar(i-1) + 24
  !       dellar(i) = 24.d0
  !       elltoi(ellar(i)) = i

  !    else
  !       ellar(i)  = ellar(i-1) + 50
  !       dellar(i) = 50.d0
  !       elltoi(ellar(i)) = i

  !    endif
  !    !write(*,*) 'ell:', ellar(i)
  ! enddo
  do i  = 1, 800
    if (i .le. 47) then
       ellar(i) = i+1
       dellar(i) = 1.d0
       elltoi(ellar(i)) = i
    elseif (i .le. 123 .and. i .ge. 48) then
       ellar(i) = ellar(i-1) + 3
       dellar(i) = 3.d0
       elltoi(ellar(i)) = i
    elseif (i .le. 190 .and. i .ge. 124) then
       ellar(i)  = ellar(i-1) + 6
       dellar(i) = 6.d0
       elltoi(ellar(i)) = i
    elseif (i .le. 450 .and. i .ge. 191) then
       ellar(i)  = ellar(i-1) + 8
       dellar(i) = 8.d0
       elltoi(ellar(i)) = i
    else
       ellar(i)  = ellar(i-1) +10
       dellar(i) = 10.d0
       elltoi(ellar(i)) = i
    endif
    !write(,) 'ell:', ellar(i)
 enddo
!!$  do i  = 1, 2048
!!$     if (i .le. 2048) then
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

!!$  do i = 1, 512
!!$     if(ellar(i) .ge. 3000) then
!!$        write(*,*) i
!!$        exit
!!$     endif
!!$  enddo 
!!$  stop
  
  lmax = 5000

  lmin = 2

  allocate(Cl(4,lmax))
  Cl(:,:) = 0d0
  allocate(Cll(4,lmax))
  Cll(:,:) = 0d0
  allocate(invCll(2,2,lmax))
  invCll(:,:,:) = 0d0
  allocate(Cllm(2,2,lmax))
  Cllm(:,:,:) = 0d0
  allocate(pClpp(3,lmax))
  pClpp(:,:) = 0d0

  Folder1 = './SOspectra/'
  Clfile = trim(Folder1)//trim('SOspectra_lenspotentialCls.dat')
  Cllfile = trim(Folder1)//trim('SOspectra_lensedCls.dat')
  !from Alex:

  !Cllfile = './SO_forecasts/CAMB/cosmo2017_10K_acc3_lensedCls.dat'



  call getenv('SCRATCHDIR',tensDir)
  !Will:
  tensDir = TRIM(tensDir)//'/Data/alphaBetaDir/'
  !Daan:
  !tensDir = '/mnt/raid-cita/meerburg/SO_forecasts/alphabeta'
  write(*,*) tensDir      
  !allocate and read bessel transforms
  !you should check the subroutine to see if your file is directed correctly
  !it is now set to my directory 
  !allocate and read bessel transforms
  !you should check the subroutine to see if your file is directed correctly
  !it is now set to my directory
  !Will
  alphabetafile = TRIM(tensDir)//'/l_r_alpha_beta.txt.MAX4000'
  alphabetaPolfile = TRIM(tensDir)//'/l_r_alpha_beta_Pol.txt.MAX4000'
  !Daan
  !alphabetafile = TRIM(tensDir)//'l_r_alpha_beta_new_Lmax5000.txt'
  !alphabetaPolfile = TRIM(tensDir)//'l_r_gamma_delta_new_Lmax5000.txt'
  !P%flagDoWigner=flagDoWigner
  if (shape .eq. 5) then
     P%flagDoWigner = 5
  endif

  P%flagDoWigner=0
  P_ISWlens%flagDoWigner = 5
  call allocate_besseltransforms(P,alphabetafile,alphabetaPolfile,cllFile,ClFile)
  call allocate_besseltransforms(P_ISWlens,alphabetafile,alphabetaPolfile,cllFile,ClFile)
  P_ISWlens%flagDoWigner = 0

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

     detCovCV = (Cll(1,j)*Cll(2,j)-Cll(4,j)**2) 
     !inverse covariance in CV limit
     if (nfields .eq. 1) then
        invCll(1,1,j) = 1.d0/Cll(1,j)
     elseif (minfields .eq. 2) then
        invCll(2,2,j) = 1.d0/Cll(2,j)
     else
        invCll(1,1,j) = Cll(2,j)/detCovCV
        invCll(1,2,j) = -Cll(4,j)/detCovCV
        invCll(2,1,j) = -Cll(4,j)/detCovCV
        invCll(2,2,j) = Cll(1,j)/detCovCV
     endif

     Cllm(1,1,j) = Cll(1,j)!/detCovCV
     Cllm(1,2,j) = Cll(4,j)!/detCovCV
     Cllm(2,1,j) = Cll(4,j)!/detCovCV
     Cllm(2,2,j) = Cll(2,j)!/detCovCV

     !write(*,*) l1,Cll(1:4,j) 
  enddo
  close(17)
  close(18)

  DB = 0.d0
  !structure
  !l2/l3/l3/l2'/l3'
  lmax = 200
  lmin = 2

  !intmax = 60 : lmax = 100
  !intmax = 210 : lmax = 3810
  intmax = 100
  imin = 40
  lmax = ellar(intmax)
  lmin = ellar(imin)

  nPhi = 200

  !to make it faster:
  stp  = 1

  ! if (nfields .eq.2) then
  !   if (shape .eq. 1) then
  !     open(unit=12,file='BispectrumCovariance_local_Pol_3.0.dat', status = 'replace')
  !   elseif(shape.eq.2) then
  !     open(unit=12,file='BispectrumCovariance_equil_Pol_3.0.dat', status = 'replace')
  !   else
  !     open(unit=12,file='BispectrumCovariance_orth_Pol_3.0.dat', status = 'replace')
  !   endif
  ! else
  !   if (shape .eq. 1) then
  !     open(unit=12,file='BispectrumCovariance_local_3.0.dat', status = 'replace')
  !   elseif(shape.eq.2) then
  !     open(unit=12,file='BispectrumCovariance_equil_3.0.dat', status = 'replace')
  !   else
  !     open(unit=12,file='BispectrumCovariance_orth_3.0.dat', status = 'replace')
  !   endif
  ! endif
  
  close(12)
  do lmaxTmp = 600, 3800,200
    do intmax = imin,size(ellar)
        if (ellar(intmax) .gt. lmaxTmp) then
          exit
      endif
    enddo

    intmax= intmax-1
    lmax = ellar(intmax)
    SumGauss =0.d0
    SumNGauss = 0.d0
    SumNGaussAll = 0.d0
    write(*,*) 'lmax:',lmax
    write(*,*) 'lmin:',lmin
    write(*,*) 'nfields:',nfields

    TotSumCV = 0.d0
    TotSumCVISWLens = 0.d0
    TotSumCVLensCross = 0.d0
    !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(dynamic) &
    !$OMP PRIVATE(l1,l2,l3, min_l,max_l,m1,p1,q1,m2,p2,q2), &
    !$OMP PRIVATE(Det,TempCovCV,atj,bis,bis_ISWlens,tmpPrefac) &
    !$OMP PRIVATE(DetISWLens,DetLensCross,fnlISW,fnl,bispectrum,bispectrum_ISWlens) &
    !$OMP REDUCTION(+:TotSumCV,TotSumCVISWLens,TotSumCVLensCross) 
    do l1 = 2, lmax
       allocate(bis(nfields**3,lmax))
       allocate(bis_ISWlens(nfields**3,lmax))
       allocate(bispectrum(nfields,nfields,nfields,1,lmax))
       allocate(bispectrum_ISWlens(nfields,nfields,nfields,1,lmax))       
       do l2 =  max(2,l1), lmax
          min_l = max(abs(l1-l2),l2)
          if (mod(l1+l2+min_l,2)/=0) then
             min_l = min_l+1 !l3 should only lead to parity even numbers
          end if
          max_l = min(lmax,l1+l2)

          bis = 0.d0
          bis_ISWlens = 0.d0
          bispectrum = 0.d0
          bispectrum_ISWlens = 0.d0
          call get_bispectrum_sss(P,l1,l2,2,lmax,shape,nfields,bis)
          call get_bispectrum_sss(P_ISWlens,l1,l2,2,lmax,5,nfields,bis_ISWlens)
          call GetThreeJs(atj(abs(l2-l1)),l1,l2,0,0)
          call reshapeBispectrum(bis,bispectrum,1,minfields,nfields)
          call reshapeBispectrum(bis_ISWlens,bispectrum_ISWlens,1,minfields,nfields)

          do l3=min_l,max_l, 2 !sum has to be even
             tmpPrefac = (atj(l3)*prefactor(l1,l2,l3)*.5)**2/tr(l1,l2,l3)
             do m2  = minfields,nfields !T,E (8 terms only)
                do p2 = minfields,nfields !T,E
                   do q2 = minfields,nfields !T,E
                      do m1  = minfields,nfields !T,E (8 terms only)
                         do p1 = minfields,nfields !T,E
                            do q1 = minfields,nfields !T,E
                               TotSumCVLensCross = TotSumCVLensCross+bispectrum_ISWlens(m1,p1,q1,1,l3)*invCll(m1,m2,l1)* &
                                    invCll(p1,p2,l2)*invCll(q1,q2,l3)*bispectrum(m2,p2,q2,1,l3)*tmpPrefac
                               TotSumCVISWLens = TotSumCVISWLens+bispectrum_ISWlens(m1,p1,q1,1,l3)*invCll(m1,m2,l1)* &
                                    invCll(p1,p2,l2)*invCll(q1,q2,l3)*bispectrum_ISWlens(m2,p2,q2,1,l3)*tmpPrefac
                               TotSumCV = TotSumCV +bispectrum(m1,p1,q1,1,l3)*invCll(m1,m2,l1)* &
                                    invCll(p1,p2,l2)*invCll(q1,q2,l3)*bispectrum(m2,p2,q2,1,l3)*tmpPrefac
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo !l3 loop
       enddo !l2 loop
       deallocate(bis)
       deallocate(bis_ISWlens)
       deallocate(bispectrum)
       deallocate(bispectrum_ISWlens)
    enddo !L1 loop
    !$OMP END PARAllEl DO


    DetFishCV = TotSumCVISWLens*TotSumCV -TotSumCVLensCross**2
    alpha = TotSumCVISWLens/DetFishCV !C/det
    beta = -TotSumCVLensCross/DetFishCV !A/det 
    if(.not. want_ISW_correction) then
       alpha=1./TotSumCV
       beta = 0.
    endif
    write(*,*) 'lensing-ISW-fnl_local correlation coefficient:', TotSumCVLensCross/TotSumCV**(1./2)/TotSumCVISWLens**(1./2)
    write(*,*) 'Fisher local error:', 1/TotSumCV**(1./2)
    write(*,*) 'Fisher ISW-lensing error:', 1/TotSumCVISWLens**(1./2)
    write(*,*) 'alpha', alpha, 'beta', beta
    !endif

    !lmax = ellar(intmax)
    do i = imin, intmax !l2 loop
       !do i = 2, 20   
       l1a = ellar(i)
       l1b = l1a

       allocate(bispectrum(nfields,nfields,nfields,lmax,lmax))
       allocate(bispectrum_ISWlens(nfields,nfields,nfields,lmax,lmax))
       allocate(bis(nfields**3,lmax))
       allocate(bis_ISWlens(nfields**3,lmax))
       bispectrum(:,:,:,:,:) = 0.d0
       bispectrum_ISWlens(:,:,:,:,:) =0.d0

       do l2 = lmin, lmax
          bis = 0.d0
          call get_bispectrum_sss(P,l1a,l2,2,lmax,shape,nfields,bis)
          ! bispectrum(1,1,1,l2,:) = bis(1,:)*.5 ! .5 for account for 2 in prefactor
          call reshapeBispectrum(bis,bispectrum,l2,minfields,nfields)
          bis_ISWlens = 0d0
          call get_bispectrum_sss(P_ISWlens,l1a,l2,2,lmax,5,nfields,bis_ISWlens)
          call reshapeBispectrum(bis_ISWlens*.5,bispectrum_ISWlens,l2,minfields,nfields)

          min_l = max(abs(l1a-l2),lmin) 
          max_l = min(lmax,l1a+l2)
          do l3=min_l,max_l,1 !sum has to be even

            call applyInvC(bispectrum,invCll,minfields,nfields, l1a,l2,l3)
            call applyInvC(bispectrum_ISWlens,invCll,minfields,nfields,l1a,l2,l3)
             !write(*,'(3I4,2E18.7)') l1,l2,l3,bispectrum(1,1,1,l2,l3),bispectrum(2,2,2,l2,l3)
          enddo
          ! write(*,*),bispectrum_ISWlens(1,1,1,l2,l3),bispectrum_ISWlens(1,1,1,l3,l2)

          ! bispectrum_ISWlens(1,l2,:) = bis_ISWlens(1,:)*.5 ! .5 for account for 2 in prefactor
       enddo

       do l2 = lmin, lmax
          max_l = min(lmax,l1a+l2)
          min_l = max(abs(l1a-l2),lmin)  
          do l3 = min_l,max_l,1
             call permuteBis(bispectrum_ISWlens,minfields,nfields,l2,l3)
             call permuteBis(bispectrum,minfields,nfields,l2,l3)
             !write(*,*) l2,l3 bispectrum(1,1,1,l2,l3),bispectrum(1,1,1,l3,l2)
             ! bispectrum(1,l3,l2) = bispectrum(1,l2,l3)
             ! bispectrum_ISWlens(1,l3,l2) = bispectrum_ISWlens(1,l2,l3)
             !bispectrum(1,l3,l2) = floc(l1,l2,l3)
             !bispectrum(1,l2,l3) = floc(l1,l2,l3)
             !write(*,*) bispectrum(1,l2,l3),bispectrum(1,l3,l2),floc(l1,l2,l3),floc(l1,l3,l2)
          enddo
       enddo



       !$OMP PARALLEL DO DEFAUlT(SHARED),SCHEDULE(dynamic) &
       !$OMP PRIVATE(l2a,l3a,l2b,l3b,phi21a,phi21b,phi21amax,phi21bmax,phi21amin,phi21bmin,j,k,m,n), &
       !$OMP PRIVATE(fnl,fnlb,fnlISW,fnlISWb,signsq,tempfac,tempfacFcM,DB_tmp), &
       !$OMP PRIVATE(p1,q1,m1,p2,q2,m2), &
       !$OMP PRIVATE(Rl2adotl3b,Rl2bdotl3a,Rl3adotl3b ,Rl2adotl2b,Rabsl2al2b ,Rabsl2al3b), &
       !$OMP PRIVATE(phi23a,phi23b,phi31a,phi31b,phi2a3b,phi2b3a,phi2a2b,phi3a3b,nPhia,nPhib,dPhia,dPhib), &
       !$OMP PRIVATE(l1dotl2a,l1adotl2b,l2dotl3a,l2dotl3b,l2adotl3b,l2bdotl3a,l3adotl3b,l2adotl2b), &
       !$OMP PRIVATE(absl2al3b,absl2al2b,absl3al3b,absl2bl3a,DB, DSNonGauss, DSNGauss, measureG, measureNG), &
       !$OMP REDUCTION(+:SumGauss,SumNGauss,SumNGaussAll)
       do j = imin, intmax !l3 loop
          !do l3a = lmin, lmax   
          l2a = ellar(j)
          !set minimum and max values of third leg (triangle constraint)
          !l3almax = min(l2a+l1a,lmax)
          !l3almin = max(abs(l2a-l1a),lmin)

          !number of steps
          nPhia = nPhi!min(l3almax-l3almin,nPhi)

          !phi21amin = -pi
          phi21amin = 0.d0
          if (abs(l2a-l1a) .gt. lmin) then
             phi21amin = 0.d0
          else
             phi21amin = angle(l1a,l2a,lmin)
          endif
          phi21amax = pi

          !deltaphi
          dPhia = (phi21amax-phi21amin)/nPhia

          do m = 0, max(nPhia,0) !phi21a loop
             phi21a = dPhia*m+phi21amin

             l3a = nint(length(l1a,l2a,phi21a))
             if ((l3a .lt. lmin) .or. (l3a .gt. lmax)) then
                !   write(*,*), lmin,lmax,l1a,l2a,l3a
                cycle
             endif
             l1dotl2a = l1a*l2a*cos(pi - phi21a)

             if(l1dotl2a .ne. l1dotl2a) write(*,*) 'NANs'

             !other angles; needed for 11 permutations 
             phi31a = angle(l3a,l1a,l2a)
             !phi31a = pi - phi23a - phi21a

             !get the SW bispectrum for the triplet (l1, l2 ,l3)
             !fnl = floc(l1a, l2a,l3a)

             !Eauate the l1 and l1'
             

             do k = imin, intmax

                l2b = ellar(k)

                nPhib = nPhi

                !phi21bmin = -pi
                phi21bmin = 0.d0
                if (abs(l2b-l1b) .gt. lmin) then
                   phi21bmin = 0
                else
                   phi21bmin = angle(l1b,l2b,lmin)
                endif
                phi21bmax = pi
                
                dPhib = (phi21bmax-phi21bmin)/nPhib
  !!$              
                do n = 0, max(nPhib,0)

                   phi21b = dPhib*n+phi21bmin
                   l3b = nint(length(l1b,l2b,phi21b))
                   if ((l3b .lt. lmin) .or. (l3b .gt. lmax)) then
                      cycle
                   endif
                   l1adotl2b = l1a*l2b*cos(pi - phi21b)
                   if(l1adotl2b .ne. l1adotl2b) write(*,*) 'NANs'
                   !prime triangle angle between l2' and l3'
                   phi31b = angle(l3b,l1b,l2b)
                   !phi31b = pi - phi23b - phi21b

                   !compute fnl^2



                   !other angles; needed for 11 permutations 
                   !phi21b = angle(l1b,l2b,l3b)                
                   if (want_all_terms) then
                     phi2a3b =  phi21a + phi31b !correct direction
                     phi2b3a =  phi21b + phi31a
                     phi2a2b = phi21a - phi21b
                     phi3a3b = phi31a - phi31b

                     !other inner products:
                     l2adotl3b = l2a*l3b*cos(phi2a3b) !l2.l3'
                     l2bdotl3a = l2b*l3a*cos(phi2b3a) !l2'.l3

                     l3adotl3b = l3a*l3b*cos(phi3a3b) !l3.l3'
                     l2adotl2b = l2a*l2b*cos(phi2a2b) !l2.l2'

                     ! !other lengths:
                     absl2al2b = nint(length(l2a,l2b,phi2a2b)) !|l2+l2'|
                     absl2al3b = nint(length(l2a,l3b,phi2a3b)) !|l3+l3'|  

                     phi2a3b =  phi21a - phi31b !correct direction
                     phi2b3a =  phi21b - phi31a
                     phi2a2b = phi21a + phi21b
                     phi3a3b = phi31a + phi31b

                     !other inner products:
                     Rl2adotl3b = l2a*l3b*cos(phi2a3b) !l2.l3'
                     Rl2bdotl3a = l2b*l3a*cos(phi2b3a) !l2'.l3

                     Rl3adotl3b = l3a*l3b*cos(phi3a3b) !l3.l3'
                     Rl2adotl2b = l2a*l2b*cos(phi2a2b) !l2.l2'

                     ! !other lengths:
                     Rabsl2al2b = nint(length(l2a,l2b,phi2a2b)) !|l2+l2'|
                     Rabsl2al3b = nint(length(l2a,l3b,phi2a3b)) !|l3+l3'|  
                   endif               

                   ! absl2al3b = nint(length(l2a,l3b,phi2a3b)) !|l2+l3'|
                   ! absl2bl3a = nint(length(l2b,l3a,phi2b3a)) !|l2'+l3|

                   !if(absl2bl3a .lt. lmin) absl2bl3a = lmin
                   !if(absl2al2b .lt. lmin) absl2al2b = lmin

                   !if(absl2bl3a .gt. lmax) absl2bl3a = lmax
                   !if(absl2al2b .gt. lmax) absl2al2b = lmax
  !!$                 
                   !3 unique terms 
                   ! DB(1) =  Cll(1,l1a)*Cll(1,l2a)*Cll(1,l2b)*pClpp(1,l1a)* &
                   !     (l1adotl2b*l1dotl2a)
                   DB(2) = 0.d0
                   DB(3) = 0.d0
                   if (want_all_terms) then
                     if(absl2al3b .gt. lmin .and. absl2al3b .lt. 4000) then                  
                         DB(2) = ( pClpp(1,nint(absl2al3b))* &
                              (l2a**2-l2adotl3b)*(l2b**2-l2bdotl3a))
                    endif
                    if(Rabsl2al3b .gt. lmin .and. Rabsl2al3b .lt. 4000) then
                         DB(2) = DB(2)+ ( pClpp(1,nint(Rabsl2al3b))* &
                              (l2a**2-Rl2adotl3b)*(l2b**2-Rl2bdotl3a))
                     endif
                     if(absl2al2b .gt. lmin .and. absl2al2b .lt. 4000) then                
                        DB(3) = (pClpp(1,nint(absl2al2b))* &
                              (l3b**2-l3adotl3b)*(l2b**2-l2adotl2b))  
                    endif
                    if(Rabsl2al2b .gt. lmin .and. Rabsl2al2b .lt. 4000) then
                        DB(3) =  DB(3)+(pClpp(1,nint(Rabsl2al2b))* &
                              (l3b**2-Rl3adotl3b)*(l2b**2-Rl2adotl2b))      
                     endif
                   endif
                   ! DB(2) = 0.d0
                   ! DB(3) = -DB(3)
                   ! write(*,*)  (l2a**2-l2adotl3b)/l2a**2,(l2b**2-l2bdotl3a)/l2b**2,(l2a**2-Rl2adotl3b)/l2a**2,(l2b**2-Rl2bdotl3a)/l2b**2,(l3b**2-l3adotl3b)/l3b**2,(l2b**2-l2adotl2b)/l2b**2,(l3b**2-Rl3adotl3b)/l3b**2,(l2b**2-Rl2adotl2b)/l2b**2
                   ! write(*,*) nint(absl2al3b),nint(absl2al2b),pClpp(1,l1a)*(l1adotl2b*l1dotl2a), pClpp(1,nint(absl2al3b))*(l2a**2-l2adotl3b)*(l2b**2-l2bdotl3a),pClpp(1,nint(Rabsl2al3b))*(l2a**2-Rl2adotl3b)*(l2b**2-Rl2bdotl3a),pClpp(1,nint(absl2al2b))*(l3b**2-l3adotl3b)*(l2b**2-l2adotl2b),(pClpp(1,nint(Rabsl2al2b))*(l3b**2-Rl3adotl3b)*(l2b**2-Rl2adotl2b)) 
                   ! 2 * pi (from phi integral) / pi - from delta 0

                   !signsq = fnl*fnlb


                   measureNG = 2*pi*dellar(i)*dellar(j)*dellar(k)*l1a*l2a*dPhia*l2b*dPhib
                   tempfacFcM(:,:) = 0
                   tempfac =  4.*measureNG*pClpp(1,l1a)*(l1adotl2b*l1dotl2a)/(2.*pi)**4/pi
                   tempfacFcM(1,1) = tempfac
                   if (nfields .gt. 1) then
                     tempfacFcM(2,1) = cos(2*(phi21a+phi31a))*tempfac
                     tempfacFcM(1,2) = cos(2*(phi21b+phi31b))*tempfac
                     tempfacFcM(2,2) = cos(2*(phi21a+phi31a))*cos(2*(phi21b+phi31b))*tempfac
                     !write(*,*) tempfacFcM(2,2)
                   endif
                   DSNonGauss = 0
                   do m2  = minfields,nfields !T,E (8 terms only)
                     do p2 = minfields,nfields !T,E
                        do q2 = minfields,nfields !T,E
                           do m1  = minfields,nfields !T,E (8 terms only)
                              do p1 = minfields,nfields !T,E
                                 do q1 = minfields,nfields !T,E
                                    DB_tmp = Cllm(m1,m2,l1a)*Cllm(p1,q1,l2a)*Cllm(p2,q2,l2b)*tempfacFcM(q1,q2)
                                    fnl = bispectrum(m1,p1,q1,l2a,l3a)
                                    fnlb = bispectrum(m2,p2,q2,l2b,l3b)
                                    fnlISW = bispectrum_ISWlens(m1,p1,q1,l2a,l3a)
                                    fnlISWb = bispectrum_ISWlens(m2,p2,q2,l2b,l3b)

                              !write(*,*) bispectrum(2,2,2,l2,l3),bispectrum(1,1,1,l2,l3)
                                    signsq = ( alpha**2*fnl*fnlb &
                                         + alpha*beta*fnl*fnlISWb + alpha*beta*fnlb*fnlISW + beta**2*fnlISW*fnlISWb)
                                    DSNonGauss =DSNonGauss+ signsq*DB_tmp
                                 enddo
                              enddo
                          enddo
                        enddo
                     enddo
                   enddo
                   

                   SumNGauss =  SumNGauss + DSNonGauss
                   if (want_all_terms) then
                       measureNG = 2*pi*dellar(i)*dellar(j)*dellar(k)*l1a*l2a*dPhia*l2b*dPhib
                       tempfacFcM(:,:) = 0
                       tempfac =  2*measureNG/(2.*pi)**4/pi
                       tempfacFcM(1,1) = tempfac
                       if (nfields .gt. 1) then
                         tempfacFcM(2,1) = cos(2*(phi21a+phi31a))*tempfac
                         tempfacFcM(1,2) = cos(2*(phi21b+phi31b))*tempfac
                         tempfacFcM(2,2) = cos(2*(phi21a+phi31a))*cos(2*(phi21b+phi31b))*tempfac
                         !write(*,*) tempfacFcM(2,2)
                       endif
                       DSNonGauss = 0
                       do m2  = minfields,nfields !T,E (8 terms only)
                         do p2 = minfields,nfields !T,E
                            do q2 = minfields,nfields !T,E
                               do m1  = minfields,nfields !T,E (8 terms only)
                                  do p1 = minfields,nfields !T,E
                                     do q1 = minfields,nfields !T,E
                                        DB_tmp = Cllm(m1,m2,l1a)*Cllm(p1,q1,l2a)*Cllm(p2,q2,l2b)*tempfacFcM(q1,q2)*DB(2)
                                        DB_tmp = DB_tmp + Cllm(m1,m2,l1a)*Cllm(p1,q1,l2a)*Cllm(p2,q2,l3a)*tempfacFcM(p2,q2)*DB(3)
                                        fnl = bispectrum(m1,p1,q1,l2a,l3a)
                                        fnlb = bispectrum(m2,p2,q2,l2b,l3b)
                                        fnlISW = bispectrum_ISWlens(m1,p1,q1,l2a,l3a)
                                        fnlISWb = bispectrum_ISWlens(m2,p2,q2,l2b,l3b)

                                  !write(*,*) bispectrum(2,2,2,l2,l3),bispectrum(1,1,1,l2,l3)
                                        signsq = ( alpha**2*fnl*fnlb &
                                             + alpha*beta*fnl*fnlISWb + alpha*beta*fnlb*fnlISW + beta**2*fnlISW*fnlISWb)
                                        DSNonGauss =DSNonGauss+ signsq*DB_tmp
                                     enddo
                                  enddo
                              enddo
                            enddo
                         enddo
                       enddo
                       SumNGaussAll =  SumNGaussAll + DSNonGauss
                   endif

                enddo !phi12b

             enddo !l2b
             DSNGauss = 0
             measureG = 2*pi*dellar(i)*dellar(j)*l1a*l2a*dPhia
             do m2  = minfields,nfields !T,E (8 terms only)
               do p2 = minfields,nfields !T,E
                  do q2 = minfields,nfields !T,E
                     do m1  = minfields,nfields !T,E (8 terms only)
                        do p1 = minfields,nfields !T,E
                           do q1 = minfields,nfields !T,E
                               fnl = bispectrum(m1,p1,q1,l2a,l3a)
                               fnlb = bispectrum(m2,p2,q2,l2a,l3a)
                               fnlISW = bispectrum_ISWlens(m1,p1,q1,l2a,l3a)
                               fnlISWb = bispectrum_ISWlens(m2,p2,q2,l2a,l3a)

                              !write(*,*) bispectrum(2,2,2,l2,l3),bispectrum(1,1,1,l2,l3)
                               signsq = ( alpha**2*fnl*fnlb &
                                   + alpha*beta*fnl*fnlISWb + alpha*beta*fnlb*fnlISW + beta**2*fnlISW*fnlISWb)
                               DSNGauss = DSNGauss + 2.*measureG*signsq*Cllm(m1,m2,l1a)*Cllm(p1,p2,l2a)*Cllm(q1,q2,l3a)/6./(2.*pi)**2/pi 
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
             enddo
             SumGauss = SumGauss + DSNGauss
             SumNGauss =  SumNGauss + DSNGauss
             SumNGaussAll =  SumNGaussAll  + DSNGauss
          enddo !phi12a
       enddo !l3a
       !$OMP END PARAllEl DO

       deallocate(bis)
       deallocate(bispectrum)
       deallocate(bis_ISWlens)
       deallocate(bispectrum_ISWlens)
       write(*,'(I4,6E16.7)') l1a,SumNGauss/alpha,SumGauss/alpha, alpha/(SumNGauss-SumGauss+alpha),SumNGaussAll/alpha,(SumNGaussAll-SumGauss)/(SumNGauss-SumGauss), alpha/(SumNGaussAll-SumGauss+alpha)


      !write(*,'(A5,I4,6E18.7)') 'l1 = ', l1a, SumGauss, SumNGauss, SumNGaussAll,TotSumCV/(TotSumCV+SumNGaussAll-SumGauss), sqrt(SumGauss/SumNGauss), sqrt(SumGauss/SumNGaussAll)
    enddo !l1a
    if (want_ISW_correction) then
      suffix = '_wISWmarg.dat'
    else
      suffix = '.dat'
    endif
    if (nfields .eq.2) then
      if (shape .eq. 1) then
        open(unit=12,file='BispectrumCovariance_local_Pol_3.0'//trim(suffix),action='write',position='append')
      elseif(shape.eq.2) then
        open(unit=12,file='BispectrumCovariance_equil_Pol_3.0'//trim(suffix),action='write',position='append')
      else
        open(unit=12,file='BispectrumCovariance_orth_Pol_3.0'//trim(suffix),action='write',position='append')
      endif
    else
      if (shape .eq. 1) then
        open(unit=12,file='BispectrumCovariance_local_3.0'//trim(suffix),action='write',position='append')
      elseif(shape.eq.2) then
        open(unit=12,file='BispectrumCovariance_equil_3.0'//trim(suffix),action='write',position='append')
      else
        open(unit=12,file='BispectrumCovariance_orth_3.0'//trim(suffix),action='write',position='append')
      endif
    endif
    write(*,'(I4,6E18.7)')   ellar(intmax),1/alpha, sqrt(SumGauss/SumNGauss), (SumNGauss-SumGauss)/TotSumCV, SumNGauss/alpha,SumNGauss/alpha, alpha/(SumNGauss-SumGauss+alpha)
    write(12,'(I4,6E18.7)') ellar(intmax), 1/alpha, sqrt(SumGauss/SumNGauss), (SumNGauss-SumGauss)/TotSumCV, SumNGauss/alpha,SumNGauss/alpha, alpha/(SumNGauss-SumGauss+alpha)
    write(*,'(A12,X,I4,X,A19,X,F11.3,X, A27,X,F11.4,A1)') 'For lmax = ',  lmax, 'the error on fnl = ', sqrt(1./alpha), 'error with full sky result:', 100*sqrt(SumGauss/alpha)-100,'%'
    ! write(*,'(I4,5E17.8)') ellar(intmax), SumGauss, SumNGauss, SumNGaussAll, sqrt(SumGauss/SumNGauss), sqrt(SumGauss/SumNGaussAll)
    !write(*,'(A12,X,I4,X,A19,X,F11.3,X, A27,X,F11.4,A1)') 'For lmax = ',  lmax, 'the error on fnl = ', sqrt(1./SumGauss), 'error with full sky result:', 100*sqrt(TotSumCV/SumGauss)-100,'%'
    !write(12,'(A5,I4,6E18.7)') 'l1 = ', l1a, SumGauss, SumNGauss, SumNGaussAll,TotSumCV/(TotSumCV+SumNGaussAll-SumGauss), sqrt(SumGauss/SumNGauss), sqrt(SumGauss/SumNGaussAll)
    close(12)
  enddo

  deallocate(Cl, Cll, pClpp,invCll,Cllm)
contains



  ! real(dl) function floc(l1,l2,l3)
  !   !SW approximation 
  !   integer :: l1, l2, l3
  !   real(dl) :: amp
  !   real(dl) :: As = 2.1056d-9
  !   !from https://arxiv.org/pdf/0812.3413.pdf Eq. 19 and 20
  !   amp = (2.d0/27./pi**2)*As
  !   floc = 1.d0/(l1+1.d0)/l1/l2/(l2+1.d0) + 1.d0/(l3+1.d0)/l3/l2/(l2+1.d0) + &
  !        1.d0/(l1+1.d0)/l1/l3/(l3+1.d0)
  !   floc = floc*amp*2.E-7 !2.E-7  is introduced to get roughly same amplitude at l_max = 500 to full fnl_local
  ! end function floc

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


  subroutine reshapeBispectrum(bis,bispectrum,l2,minfields,nfields)
    real(dl), intent(in) :: bis(:,:)
    integer, intent(in) :: nfields,minfields
    real(dl) :: bispectrum(:,:,:,:,:)
    integer ::  l2
    if (nfields .eq. 1) then
       bispectrum(1,1,1,l2,:) = bis(1,:)
    elseif (minfields .eq. 2) then
       bispectrum(2,2,2,l2,:) = bis(8,:)
    else
       bispectrum(1,1,1,l2,:) = bis(1,:)
       bispectrum(1,1,2,l2,:) = bis(2,:)
       bispectrum(1,2,1,l2,:) = bis(3,:)
       bispectrum(1,2,2,l2,:) = bis(4,:)
       bispectrum(2,1,1,l2,:) = bis(5,:)
       bispectrum(2,1,2,l2,:) = bis(6,:)
       bispectrum(2,2,1,l2,:) = bis(7,:)
       bispectrum(2,2,2,l2,:) = bis(8,:)
    endif
  end subroutine reshapeBispectrum

  subroutine permuteBis(bispectrum,minfields,nfields,l2,l3)
    real(dl)  bispectrum(:,:,:,:,:)
    integer l2,l3,nfields,minfields
    if (nfields .eq. 1) then
       bispectrum(1,1,1,l3,l2) = bispectrum(1,1,1,l2,l3)
    elseif (minfields .eq. 2) then
       bispectrum(2,2,2,l3,l2) = bispectrum(2,2,2,l2,l3)
    else
       bispectrum(1,1,1,l3,l2) = bispectrum(1,1,1,l2,l3)
       bispectrum(1,2,1,l3,l2) = bispectrum(1,1,2,l2,l3)
       bispectrum(1,1,2,l3,l2) = bispectrum(1,2,1,l2,l3)
       bispectrum(2,1,1,l3,l2) = bispectrum(2,1,1,l2,l3)
       bispectrum(2,1,2,l3,l2) = bispectrum(2,2,1,l2,l3)
       bispectrum(2,2,1,l3,l2) = bispectrum(2,1,2,l2,l3)
       bispectrum(1,2,2,l3,l2) = bispectrum(1,2,2,l2,l3)
       bispectrum(2,2,2,l3,l2) = bispectrum(2,2,2,l2,l3)
    endif
  end subroutine permuteBis

  subroutine applyInvC(bispectrum,invCl,minfields,nfields,l1,l2,l3)
    real(dl), intent(inout) :: bispectrum(:,:,:,:,:)
    real(dl), intent(in):: invCl(:,:,:)
    integer, intent(in):: nfields,minfields
    integer, intent(in) :: l1,l2,l3
    integer :: m1,m2,p1,p2,q1,q2
    real(dl) :: tmpBis(nfields,nfields,nfields)
    tmpBis(:,:,:) = 0d0
    do m1 = minfields,nfields
       do m2 = minfields,nfields
          do p1 = minfields,nfields
             do p2 = minfields,nfields
                do q1 = minfields,nfields
                   do q2 = minfields,nfields
                      tmpBis(q1,p1,m1) = tmpBis(q1,p1,m1) +invCl(q1,q2,l1)*invCl(p1,p2,l2)*invCl(m1,m2,l3)*bispectrum(q2,p2,m2,l2,l3)
                  enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    bispectrum(:,:,:,l2,l3) = tmpBis(:,:,:)
  end subroutine applyInvC

  subroutine allocate_besseltransforms(P,alphabetafile,alphabetaPolfile,Cllfile,pClfile)
    character(120),intent(in) :: alphabetafile, alphabetaPolfile,Cllfile,pClfile
    type(bfs) :: P
    character(120) :: alphaTfile, betaTfile, alphaEfile, betaEfile, rfile
    character(120) :: gammaTfile, deltaTfile, gammaEfile, deltaEfile, header
    integer :: FileUnit
    character(6) :: form = '(I3)'
    integer :: ell, ellr,j
    real(dl) :: rell !for cmbs4 noise file...
    character(6) :: clmax
    !limiting r-integration range
    integer :: dlmax
    integer :: i
    real(dl)  :: tmpCll(4,2:5000)
    real(dl) :: CMB2COBEnorm = 7428350250000.d0
    real(dl), parameter :: pi = 3.14159265359
    !real(dl) :: rmin, rmax
    !integer  :: irmin, irmax, ch1, ch2
    integer  :: dimr
    dlmax = 4000 !==> File ell max (allocation)

    FileUnit  = 10
    if(dlmax .ge. 1000) form = '(I4)'
    write(clmax, form)  dlmax
    if (P%flagDoWigner .eq. 3) then

       allocate(P%Cll(4,5000))
       P%Cll(:,:) = 0d0
       ! Cllfile = trim(Folder1)//trim('SOspectra_lensedCls.dat')
       open(unit=18,file = Cllfile, status='old')
       do j = 1, 5000
          !#    L    TT             EE             BB             TE 
          if (j .eq. 1) then
             read(18,*)
             cycle
          endif
          read(18,*) ell, tmpCll(1:4,j)
          P%Cll(1,j) = 2.*pi*tmpCll(1,j)/(real(ell,dl)*(real(ell,dl)+1.))/CMB2COBEnorm
          P%Cll(2,j) = 2.*pi*tmpCll(4,j)/(real(ell,dl)*(real(ell,dl)+1.))/CMB2COBEnorm
          P%Cll(3,j) = 2.*pi*tmpCll(2,j)/(real(ell,dl)*(real(ell,dl)+1.))/CMB2COBEnorm
          P%Cll(4,j) = 2.*pi*tmpCll(3,j)/(real(ell,dl)*(real(ell,dl)+1.))/CMB2COBEnorm
          !write(*,*) l1,Cll(1:4,j) 
       enddo
       close(18)
    endif
    if (P%flagDoWigner .eq. 5) then

       allocate(P%Cll(4,5000))
       P%Cll(:,:) = 0d0
       allocate(P%pCll(4,5000))
       P%pCll(:,:) = 0d0
       ! Cllfile = trim(Folder1)//trim('SOspectra_lensedCls.dat')
       open(unit=17,file = pClfile, status='old')
       open(unit=18,file = Cllfile, status='old')
       do j = 1, 5000
          !#    L    TT             EE             BB             TE             PP             TP             EP
          !#    L    TT             EE             BB             TE 
          if (j .eq. 1) then
             read(17,*)
             read(18,*)
             cycle
          endif
          read(17,*) ell, tmpCll(1:4,j),P%pCll(1:3,j)
          read(18,*) ell, tmpCll(1:4,j)
          P%pCll(1,j) = 2.*pi*P%pCll(1,j)/(real(ell,dl)*(real(ell,dl)+1.))**2
          P%pCll(2:3,j) = 2.*pi*P%pCll(2:3,j)/(real(ell,dl)*(real(ell,dl)+1.))**(3.d0/2.d0)/CMB2COBEnorm**(1./2)

          P%Cll(2,j) = 2.*pi*tmpCll(1,j)/(real(ell,dl)*(real(ell,dl)+1.))/CMB2COBEnorm
          P%Cll(3,j) = 2.*pi*tmpCll(4,j)/(real(ell,dl)*(real(ell,dl)+1.))/CMB2COBEnorm
          P%Cll(4,j) = 2.*pi*tmpCll(2,j)/(real(ell,dl)*(real(ell,dl)+1.))/CMB2COBEnorm
          !write(*,*) l1,Cll(1:4,j) 
       enddo
       close(17)
       close(18)
    endif
    !locate alpha,beta,gamma,delta files
    !alphabetafile = 'alphabeta/l_r_alpha_beta_new_Lmax'//trim(clmax)//'.txt' 
    ! gammadeltafile = 'alphabeta/l_r_gamma_delta_new_Lmax'//trim(clmax)//'.txt' 
    ! alphabetafile = '/Volumes/wcoulton/nonGaus/Data/alphaBetaDir/l_r_alpha_beta.txt.MAX4000'
    ! alphabetaPolfile = '/Volumes/wcoulton/nonGaus/Data/alphaBetaDir/l_r_alpha_beta_Pol.txt.MAX4000'
    ! alphabetafile = '/tigress/wcoulton/nonGaus/Data/alphaBetaDir/l_r_alpha_beta.txt.MAX4000'
    ! alphabetaPolfile = '/tigress/wcoulton/nonGaus/Data/alphaBetaDir/l_r_alpha_beta_Pol.txt.MAX4000'

    !alphabetafile = '/Volumes/wCoulton2TB/nonGaus/Data/alphaBetaDir/l_r_alpha_beta.txt.MAX4000'
    !alphabetaPolfile = '/Volumes/wCoulton2TB/nonGaus/Data/alphaBetaDir/l_r_alpha_beta_Pol.txt.MAX4000'

    if ((P%flagDoWigner .eq. 0) .or. (P%flagDoWigner .eq. 1) ) then
       open(unit=FileUnit,file = trim(alphabetafile), status='old')
       open(unit=FileUnit+1,file = trim(alphabetaPolfile), status='old')
       !reading in header (nice illustration of the limitations of Fortran)
       !blank line
       write(*,*)
       do i=1,3
          !alpha/beta file, show header:
          read(FileUnit,'(1A70)') header
          if(i .eq. 1 .or. i .eq. 2) write(*,*) trim(header)
       enddo
       write(*,*)
       dimr = 438

       do i=1,3
          !gamma/delta file, skip header:
          read(FileUnit+1,'(1A70)') header
       enddo

       allocate(P%rarray(dimr), P%deltar(dimr))

       allocate(P%ar(dimr,dlmax,2))
       allocate(P%br(dimr,dlmax,2))
       allocate(P%gr(dimr,dlmax,2))
       allocate(P%dr(dimr,dlmax,2))
       P%ar(:,:,:) = 0d0
       P%br(:,:,:) = 0d0
       P%gr(:,:,:) = 0d0
       P%dr(:,:,:) = 0d0
       !read
       write(*,*) 'reading alpha/beta/gamma/delta files ...'
       write(*,*) ' Will altered this. I get from camb only 438 points instead of 603 points'
       write(*,*) ' commented out some statements in get_bispectrum. Maybe put them here for clarity?'
       write(*,*) ' altered get_bispectrum now prefactor and wigners are squared.'
       do ell = 2, dlmax
          do i  = 1, dimr
             read(FileUnit,*) ellr, P%rarray(i), P%deltar(i), P%ar(i,ell,1), P%br(i,ell,1) , P%gr(i,ell,1), P%dr(i,ell,1)
             read(FileUnit+1,*) ellr, P%rarray(i), P%deltar(i), P%ar(i,ell,2), P%br(i,ell,2) , P%gr(i,ell,2), P%dr(i,ell,2)
             !polarization alpha,beta,gamma,delta need factor ((l-2)!/(l+2)!^)1/2
             !P%br(i,ell,2) = sqrt((ell-1d0)*ell*(ell+1d0)*(ell+2d0))*P%br(i,ell,2)
             !P%ar(i,ell,2) = sqrt((ell-1d0)*ell*(ell+1d0)*(ell+2d0))*P%ar(i,ell,2)
             !P%gr(i,ell,2) = sqrt((ell-1d0)*ell*(ell+1d0)*(ell+2d0))*P%gr(i,ell,2)
             !P%dr(i,ell,2) = sqrt((ell-1d0)*ell*(ell+1d0)*(ell+2d0))*P%dr(i,ell,2)
             !write(*,*) ell, rarray(i), deltar(i), ar(i,l1,1), br(i,l1,1) , ar(i,l1,2), br(i,l1,2)
          enddo
       enddo
       close(FileUnit)
       close(FileUnit+1)
       write(*,*) 'reading complete'
    endif

  end subroutine allocate_besseltransforms

  subroutine deallocate_besseltransforms(P)
    type(bfs) :: P
    if ((P%flagDoWigner .eq. 0) .or. (P%flagDoWigner .eq. 1)) then 
       deallocate(P%ar,P%br,P%gr,P%dr)
       deallocate(P%rarray, P%deltar)
    elseif (P%flagDoWigner.eq.3) then
       deallocate(P%Cll)
    elseif (P%flagDoWigner.eq.5) then
       deallocate(P%Cll)
       deallocate(P%pCll)
    endif
  end subroutine deallocate_besseltransforms

  subroutine get_bispectrum_sss(PQ,l1,l2,ellmin,ellmax,shape,nfields,bis)  
    !opening file params 

    !real(dl), allocatable:: rarray(:), deltar(:)
    integer :: dimr
    integer :: fcnt
    type(bfs),intent(in) :: PQ
    real(dl), intent(out) :: bis(:,:)
    integer,intent(in) :: nfields
    integer :: minfields 

    !bispectra
    real(dl) :: Btemp
    real(dl) :: abint
    integer :: local = 1, equil = 2, ortho = 3, folded = 4, ISWlens = 5, localISWapprox = 6!not yet implemented
    character(7) :: nshape
    integer,intent(in) :: shape

    !loop params
    integer,intent(in) :: l1, l2
    integer :: L3
    integer :: i, p, j, k, q, r, s

    integer :: ell, ellr
    real(dl) :: rell !for cmbs4 noise file...
    character(6) :: clmax
    integer, intent(in) :: ellmax,ellmin
    integer :: lmin, max_l, min_l
    integer :: FileUnit

    !limiting r-integration range
    real(dl) :: rmin, rmax
    integer :: irmin, irmax, ch1, ch2
    logical :: lim_r_samp = .false.


    !lmin is not working yet (keep 2)
    lmin = ellmin

    !write(*,'(A6,X,I4)') 'lmax:', ellmax
    !write(*,'(A6,X,I4)') 'lmin:', ellmin    
    !once tested, you can then limit integration range
    !putting new integration range below:
    !you can set this to falls, it should still work but bea little slow. 
    lim_r_samp = .true.

    minfields  = 1 !2 for only E Can also change this. Should work. 
    if (nfields .eq. 1) minfields  = 1

    ! if (nfields .eq. 1) then
    !    write(*,*) "temperature only <TTT>"
    ! elseif(minfields .eq. 2 .and. nfields .eq.2) then
    !    write(*,*) "polarization only <EEE>"
    ! else
    !    write(*,*) "temperature and polarization <TTT>, <TTE>, <TEE>, <EEE>"
    ! endif
    !tested these below. Orthogonal becomes better constrained for narrower range of r. 
    if (shape .eq. 1) then
       nshape = 'local'
       !tested using  test_rintegral = .true.
       !will make some plots/currently tested for lmax  = 2000.
       !gives 1% accuracy 
       rmin = 13805.1
       rmax = 14067.6
    elseif(shape .eq. 2) then
       nshape = 'equil'
       rmin = 13728.1
       rmax = 14137.6
    elseif(shape .eq. 3) then
       nshape = 'ortho'
       rmin = 13675.6
       rmax = 14172.6
    elseif (shape .eq. 4) then
       nshape  = 'folded'

    elseif (shape .eq. 5) then
       nshape  = 'ISW-lens'
    elseif (shape .eq. 6) then
       nshape  = 'local-ISWapprox'
    endif
    dimr = 438
    ch1 = 1
    ch2 = 1
    !if not limiting integration range:
    irmin  = dimr
    irmax = 1
    !if limiting, change irmin and irmax
    if (((PQ%flagDoWigner .eq. 0 ).or. (PQ%flagDoWigner .eq. 1)) .and. shape .ne. 5) then
       if(lim_r_samp) then
          !write(*,*) 'limiting integration in r. 99% accuracy'
          do i  = 1, dimr
             !locate rmin and rmax
             if (PQ%rarray(i) .le. rmin .and. ch1 .eq. 1) then
                irmin = i
                ch1 = 0
             endif
             if (PQ%rarray(i) .le. rmax .and. ch2 .eq. 1) then
                irmax = i
                ch2 = 0
             endif
             !write(*,*) 'using only', irmin - irmax, 'instead of', dimr, 'points' 
          enddo
          !write(*,*) 'using only', irmin - irmax, 'instead of', dimr, 'points'
       endif
    endif
    min_l = max(abs(l1-l2),l2)

    max_l = min(ellmax,l1+l2)
    !you can set the wigners to 1 here if you want them to be excluded (I think you still need 0.5 below though)



    !$OMP parallel do private(p,i, j, k,abint,Btemp,l3,fcnt) shared(PQ,irmax,irmin,bis,nfields,shape,local,equil,ortho,l1,l2,min_l,max_l)
    do l3=min_l,max_l, 1 !sum has to be even
       fcnt = 1 !counter for TTT, TEE, TTE etc
       do p  = minfields,nfields !T,E (8 terms only)
          do j = minfields,nfields !T,E
             do k = minfields,nfields !T,E
                Btemp = 0.d0
                if ((PQ%flagDoWigner .eq. 0 ).or. (PQ%flagDoWigner .eq. 1)) then
                   if (shape .eq. 5) then
                      Btemp = fPhiISW(l1,l2,l3,PQ%pCll(1+p,l1),PQ%Cll(k+j,l3),j)+ & 
                           fPhiISW(l3,l2,l1,PQ%pCll(1+k,l3),PQ%Cll(p+j,l1),j)+ & 
                           fPhiISW(l2,l1,l3,PQ%pCll(1+j,l2),PQ%Cll(k+p,l3),p)+ & 
                           fPhiISW(l3,l1,l2,PQ%pCll(1+k,l3),PQ%Cll(j+p,l2),p)+ & 
                           fPhiISW(l2,l3,l1,PQ%pCll(1+j,l2),PQ%Cll(p+k,l1),k)+ & 
                           fPhiISW(l1,l3,l2,PQ%pCll(1+p,l1),PQ%Cll(j+k,l2),k)
                      Btemp = Btemp/2.
                   else if (shape .eq. 6) then
                      Btemp = floc(l1,l2,l3)
                   else
                      do i = irmax, irmin !rmini, rmaxi !
                         if (shape .eq. local) then 
                            abint = PQ%br(i,l1,p)*PQ%br(i,l2,j)*PQ%ar(i,l3,k) + &
                                 PQ%br(i,l3,k)*PQ%br(i,l1,p)*PQ%ar(i,l2,j) + &
                                 PQ%br(i,l2,j)*PQ%br(i,l3,k)*PQ%ar(i,l1,p)
                         elseif (shape .eq. equil) then
                            abint = -3.*( PQ%br(i,l1,p)*PQ%br(i,l2,j)*PQ%ar(i,l3,k) + &
                                 PQ%br(i,l3,k)*PQ%br(i,l1,p)*PQ%ar(i,l2,j) + &
                                 PQ%br(i,l2,j)*PQ%br(i,l3,k)*PQ%ar(i,l1,p)) + &
                                 (-6.)*(PQ%dr(i,l1,p)*PQ%dr(i,l2,j)*PQ%dr(i,l3,k)) + &
                                 (+3.)*(PQ%gr(i,l1,p)*PQ%dr(i,l2,j) + PQ%gr(i,l2,j)*PQ%dr(i,l1,p))*PQ%br(i,l3,k) + &
                                 (+3.)*(PQ%dr(i,l1,p)*PQ%br(i,l2,j) + PQ%dr(i,l2,j)*PQ%br(i,l1,p))*PQ%gr(i,l3,k) + &
                                 (+3.)*(PQ%gr(i,l1,p)*PQ%br(i,l2,j) + PQ%gr(i,l2,j)*PQ%br(i,l1,p))*PQ%dr(i,l3,k) 
                         elseif (shape .eq. ortho) then
                            abint = -9.*( PQ%br(i,l1,p)*PQ%br(i,l2,j)*PQ%ar(i,l3,k) + &
                                 PQ%br(i,l3,k)*PQ%br(i,l1,p)*PQ%ar(i,l2,j) + &
                                 PQ%br(i,l2,j)*PQ%br(i,l3,k)*PQ%ar(i,l1,p)) + &
                                 (-24.)*(PQ%dr(i,l1,p)*PQ%dr(i,l2,j)*PQ%dr(i,l3,k)) + &
                                 (+9.)*(PQ%gr(i,l1,p)*PQ%dr(i,l2,j) + PQ%gr(i,l2,j)*PQ%dr(i,l1,p))*PQ%br(i,l3,k) + &
                                 (+9.)*(PQ%dr(i,l1,p)*PQ%br(i,l2,j) + PQ%dr(i,l2,j)*PQ%br(i,l1,p))*PQ%gr(i,l3,k) + &
                                 (+9.)*(PQ%gr(i,l1,p)*PQ%br(i,l2,j) + PQ%gr(i,l2,j)*PQ%br(i,l1,p))*PQ%dr(i,l3,k) 
                         elseif (shape .eq. folded) then
                         endif

                         !1D heat map of integral to determine r values to include
                         !cum_rint(i) = cum_rint(i) + rarray(i)**2*deltar(i)*abint

                         Btemp = Btemp + PQ%rarray(i)**2*PQ%deltar(i)*abint

                      enddo !r loop
                   endif
                endif
                !means you have to recompute bispectrum elements every time, but requires low memory 
                ! Need prefactor*a3j squared as that is what the binned function computes
                ! Need factor of 0.5 as prefactor includes a factor of two. 
                ! The prefactor factor of two is needed for local, equil and orth. But not ISW lensing or others!!!! 
                ! I include a factor of a half above for the ISW for simplicity
                ! One factor is needed to replace the factor missing above however the squaring introduces an extra one.

                bis(fcnt,l3)  = 2.0*Btemp

                !write(*,*) l3, bis(fcnt,l3)
                !if you write in this loop, it severaly slows down the inner loop.
                !write(*,*) BredTemp(k,j,p)
                fcnt = fcnt + 1

             enddo !T, E loop
          enddo  !T, E loop
       enddo  !T, E loop

    enddo !l3 loop
    !$OMP END parallel do




  end subroutine get_bispectrum_sss


  real(dl) function floc(l1,l2,l3)
    !SW approximation
    integer :: l1, l2, l3
    real(dl) :: amp
    real(dl), parameter :: pi = 3.14159265359
    real(dl) :: As = 2.1056d-9
    !from https://arxiv.org/pdf/0812.3413.pdf Eq. 19 and 20
    amp = (2.d0/27./pi**2)*As
    floc = 1.d0/(l1+1.d0)/l1/l2/(l2+1.d0) + 1.d0/(l3+1.d0)/l3/l2/(l2+1.d0) + &
         1.d0/(l1+1.d0)/l1/l3/(l3+1.d0)
    floc = floc*amp*2.E-7 
  end function floc


  real function tr(l1,l2,l3)
    integer, intent(in) :: l1,l2,l3
    if ((l1.eq.l2).and.(l2.eq.l3)) then
       tr  = 6.d0
    elseif ((l1.eq.l2).or.(l2.eq.l3).or.(l3.eq.l1)) then 
       tr = 2d0
    else
       tr = 1.d0
    endif

  end function tr

  ! From https://arxiv.org/pdf/1509.08107.pdf
  real(dl) function fPhiISW(l1,l2,l3,CPT_l1,CTT_l3,p2)
    integer :: l1, l2, l3,p2
    real(dl) :: CPT_l1,CTT_l3
    fPhiISW = .5*((l1+1)*l1-l2*(l2+1)+l3*(l3+1))*CPT_l1*CTT_l3
    if (p2 .eq.2) then
       fPhiISW = fPhiISW*((l1+1)*l1-l2*(l2+1)-l3*(l3+1))*((l1+1)*l1-l2*(l2+1)-l3*(l3+1)+2) +&
            (-2.0)*fPhiISW*l2*(l2+1)*l3*(l3+1)
       fPhiISW = fPhiISW / sqrt(4.*(l2-1)*l2*(l2+1)*(l2+2)*(l3-1)*l3*(l3+1)*(l3+2) )
    endif
  end function fPhiISW

  real function prefactor(l1,l2,l3)
    real(dl), parameter :: pi = 3.14159265359
    integer, intent(in) :: l1,l2,l3

    prefactor = 2.0*sqrt((1./4.)*((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.))/pi)
  end function prefactor


  subroutine GetThreeJs(thrcof,l2in,l3in,m2in,m3in)
    !Recursive evaluation of 3j symbols. Does minimal error checking on input
    !parameters.
    implicit none
    integer, parameter :: dl = KIND(1.d0)
    integer, intent(in) :: l2in,l3in, m2in,m3in
    real(dl), dimension(*) :: thrcof
    INTEGER, PARAMETER :: i8 = selected_int_kind(18)
    integer(i8) :: l2,l3,m2,m3
    integer(i8) :: l1, m1, l1min,l1max, lmatch, nfin, a1, a2

    real(dl) :: newfac, oldfac, sumfor, c1,c2,c1old, dv, denom, x, sum1,sumuni
    real(dl) :: x1,x2,x3, y,y1,y2,y3,sum2,sumbac, ratio,cnorm, sign1, thresh
    integer i,ier, index, nlim, sign2
    integer nfinp1,nfinp2,nfinp3, lstep, nstep2,n
    real(dl), parameter :: zero = 0._dl, one = 1._dl
    real(dl), parameter ::  tiny = 1.0d-30, srtiny=1.0d-15, huge = 1.d30,srhuge = 1.d15

    ! routine to generate set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3)

    ! by recursion from l1min = max(abs(l2-l3),abs(m1)) 
    !                to l1max = l2+l3
    ! the resulting 3j-coeffs are stored as thrcof(l1-l1min+1)

    ! to achieve the numerical stability, the recursion will proceed
    ! simultaneously forwards and backwards, starting from l1min and l1max
    ! respectively.
    !
    ! lmatch is the l1-value at which forward and backward recursion are
    ! matched.
    !
    ! ndim is the length of the array thrcof
    !
    ! ier = -1 for all 3j vanish(l2-abs(m2)<0, l3-abs(m3)<0 or not integer)
    ! ier = -2 if possible 3j's exceed ndim
    ! ier >= 0 otherwise

    l2=l2in
    l3=l3in
    m2=m2in
    m3=m3in
    newfac = 0
    lmatch = 0
    m1 = -(m2+m3)

    ! check relative magnitude of l and m values
    ier = 0

    if (l2 < abs(m2) .or. l3 < m3) then
       ier = -1
       ! call MpiStop('error ier = -1')
       print*, 'error ier = -1',l2,abs(m2),l3,m3
       stop
       return
    end if

    ! limits for l1
    l1min = max(abs(l2-l3),abs(m1))
    l1max = l2+l3

    if (l1min >= l1max) then
       if (l1min/=l1max) then
          ier = -1

          !call MpiStop('error ier = -1')
          print*, 'error ier = -1',l1min,l1max 
          stop
          return
       end if

       ! reached if l1 can take only one value, i.e.l1min=l1max
       thrcof(1) = (-1)**abs(l2+m2-l3+m3)/sqrt(real(l1min+l2+l3+1,dl))
       return

    end if

    nfin = l1max-l1min+1

    ! starting forward recursion from l1min taking nstep1 steps
    l1 = l1min
    thrcof(1) = srtiny
    sum1 = (2*l1 + 1)*tiny

    lstep = 1

30  lstep = lstep+1
    l1 = l1+1

    oldfac = newfac
    a1 = (l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)
    a2 = (l1+m1)*(l1-m1)*(-l1+l2+l3+1)
    newfac = sqrt(a2*real(a1,dl))
    if (l1 == 1) then
       !IF L1 = 1  (L1-1) HAS TO BE FACTORED OUT OF DV, HENCE
       c1 = -(2*l1-1)*l1*(m3-m2)/newfac
    else

       dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 + l1*(l1-1)*(m3-m2)
       denom = (l1-1)*newfac

       if (lstep > 2) c1old = abs(c1)
       c1 = -(2*l1-1)*dv/denom

    end if

    if (lstep<= 2) then

       ! if l1=l1min+1 the third term in the recursion eqn vanishes, hence
       x = srtiny*c1
       thrcof(2) = x
       sum1 = sum1+tiny*(2*l1+1)*c1*c1
       if(lstep==nfin) then
          sumuni=sum1
          go to 230
       end if
       goto 30

    end if

    c2 = -l1*oldfac/denom

    ! recursion to the next 3j-coeff x  
    x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
    thrcof(lstep) = x
    sumfor = sum1
    sum1 = sum1 + (2*l1+1)*x*x
    if (lstep/=nfin) then

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(x) >= srhuge) then

          ! REACHED IF LAST 3J-COEFFICIENT LARGER THAN SRHUGE
          ! SO THAT THE RECURSION SERIES THRCOF(1), ... , THRCOF(LSTEP)
          ! HAS TO BE RESCALED TO PREVENT OVERFLOW

          ier = ier+1
          do i = 1, lstep
             if (abs(thrcof(i)) < srtiny) thrcof(i)= zero
             thrcof(i) = thrcof(i)/srhuge
          end do

          sum1 = sum1/huge
          sumfor = sumfor/huge
          x = x/srhuge

       end if

       ! as long as abs(c1) is decreasing, the recursion proceeds towards
       ! increasing
       ! 3j-valuse and so is numerically stable. Once an increase of abs(c1) is 
       ! detected, the recursion direction is reversed.

       if (c1old > abs(c1)) goto 30

    end if !lstep/=nfin

    ! keep three 3j-coeffs around lmatch for comparison with backward recursion

    lmatch = l1-1
    x1 = x
    x2 = thrcof(lstep-1)
    x3 = thrcof(lstep-2)
    nstep2 = nfin-lstep+3

    ! --------------------------------------------------------------------------
    !
    ! starting backward recursion from l1max taking nstep2 stpes, so that
    ! forward and backward recursion overlap at 3 points 
    ! l1 = lmatch-1, lmatch, lmatch+1

    nfinp1 = nfin+1
    nfinp2 = nfin+2
    nfinp3 = nfin+3
    l1 = l1max
    thrcof(nfin) = srtiny
    sum2 = tiny*(2*l1+1)

    l1 = l1+2
    lstep=1

    do
       lstep = lstep + 1
       l1= l1-1

       oldfac = newfac
       a1 = (l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)
       a2 = (l1+m1-1)*(l1-m1-1)*(-l1+l2+l3+2)
       newfac = sqrt(a1*real(a2,dl))

       dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 +l1*(l1-1)*(m3-m2)

       denom = l1*newfac
       c1 = -(2*l1-1)*dv/denom
       if (lstep <= 2) then

          ! if l2=l2max+1, the third term in the recursion vanishes

          y = srtiny*c1
          thrcof(nfin-1) = y
          sumbac = sum2
          sum2 = sum2 + tiny*(2*l1-3)*c1*c1

          cycle

       end if

       c2 = -(l1-1)*oldfac/denom

       ! recursion to the next 3j-coeff y
       y = c1*thrcof(nfinp2-lstep)+c2*thrcof(nfinp3-lstep)

       if (lstep==nstep2) exit

       thrcof(nfinp1-lstep) = y
       sumbac = sum2
       sum2 = sum2+(2*l1-3)*y*y

       ! see if last unnormalised 3j-coeff exceeds srhuge
       if (abs(y) >= srhuge) then

          ! reached if 3j-coeff larger than srhuge so that the recursion series
          ! thrcof(nfin),..., thrcof(nfin-lstep+1) has to be rescaled to prevent
          ! overflow

          ier=ier+1
          do i = 1, lstep
             index=nfin-i+1
             if (abs(thrcof(index)) < srtiny) thrcof(index)=zero
             thrcof(index) = thrcof(index)/srhuge
          end do

          sum2=sum2/huge
          sumbac=sumbac/huge

       end if

    end do

    ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the 
    ! corresponding backward recursion vals y1, y2, y3

    y3 = y
    y2 = thrcof(nfinp2-lstep)
    y1 = thrcof(nfinp3-lstep)

    ! determine now ratio such that yi=ratio*xi (i=1,2,3) holds with minimal
    ! error

    ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
    nlim = nfin-nstep2+1

    if (abs(ratio) >= 1) then

       thrcof(1:nlim) = ratio*thrcof(1:nlim) 
       sumuni = ratio*ratio*sumfor + sumbac

    else

       nlim = nlim+1
       ratio = 1/ratio
       do n = nlim, nfin
          thrcof(n) = ratio*thrcof(n)
       end do
       sumuni = sumfor + ratio*ratio*sumbac

    end if
    ! normalise 3j-coeffs

230 cnorm = 1/sqrt(sumuni)

    ! sign convention for last 3j-coeff determines overall phase

    sign1 = sign(one,thrcof(nfin))
    sign2 = (-1)**(abs(l2+m2-l3+m3))
    if (sign1*sign2 <= 0) then
       cnorm = -cnorm
    end if
    if (abs(cnorm) >= one) then
       thrcof(1:nfin) = cnorm*thrcof(1:nfin)
       return
    end if

    thresh = tiny/abs(cnorm)

    do n = 1, nfin
       if (abs(thrcof(n)) < thresh) thrcof(n) = zero
       thrcof(n) = cnorm*thrcof(n)
    end do
    return 

  end subroutine GetThreeJs

end program FlatSky
