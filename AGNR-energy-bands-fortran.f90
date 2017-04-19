c This fortran program uses the tight binding method to obtain the energy eigen  values for 
c arm chair graphene nanoribbon with 12 atoms in the unit cell.



    
      program ngnr


      real a,a1,a2,a3
      real pi, e2p,t,incremnt
      real, allocatable :: kx(:),ky(:)
      
      integer,allocatable :: fnn(:,:)
      real,allocatable :: fnnxc(:,:), fnnyc(:,:),x(:,:),y(:,:)
      
      
      integer i,j,kl,index1,xx,yy,fnnj,snnj,tnnj
      real g1,g2,g3,dkx,dky
      real xcatm,ycatm,xdist,ydist,distmag
      integer N,L,mind,nind,fnnatm,snnatm,tnnatm 
      integer kkx,kky
      real denomfnn,numfnn

      
      real lowlimitx,uplimitx,vrangex,incremntx
      real lowlimity,uplimity,vrangey,incremnty 
      integer nkx,nky,nkz,numpoints, count1,natoms
      real ehigh, elow, del,dinc
      real kxc,kyc, count2, eps

c     LAPACK related 
c     complex*16, dimension(10,10) :: H,S
      complex*16, allocatable :: H(:,:),S(:,:),b(:),WORK(:)
      complex*16 DUMMY(1,1)
      integer ii, ok

      open (unit = 252, file = "ngnreig.txt")
      open (unit = 253, file = "points.txt")
  
      pi = 3.14159265359
      a = 2.46 
      a1 = 0.71
      a2 = 1.2298
      a3 = 1.42
      g1 = 1.0


c Half the number of atoms in the unit cell
      N = 7
      L = 3
      natoms = 2*N
    
      allocate(x(L,natoms))
      allocate(y(L,natoms))
      allocate(H(natoms,natoms))
      allocate(S(natoms,natoms))
      allocate(b(natoms))
      allocate(WORK(2*natoms))
    
    
    
    do i = 1,N
    
    if (mod(i,2) .eq. 1) then
        x(1,i) = .71
          else
        x(1,i) = 0
          endif
    
          y(1,i) = (i-1)*1.23
          enddo


      index1 = 1
c creating the second line (mirror image of the first one)
c this completes the agnr unit cell for agnr
      do i = N +1 , natoms
    
      if (mod(i,2) .eq. 0) then
      x(1,i) = 2.13
      else
      x(1,i)= 2.84
      end if
    
      y(1,i) = y(1,i - index1)
    
      index1 = index1 + 2
   
      end do

c repeating the unit cell L times to complete 
c the structure
   
      do i = 2,L
      do j = 1, natoms
      x(i,j) = x(1,j) + 4.26*(i-1)
      y(i,j) = y(1,j)
      enddo
      enddo

      do i = 1,L
      do j = 1, natoms
      write(253,*) x(i,j), y(i,j)     
    
      enddo
      enddo

      
      eps = 0.01
      allocate(fnn(natoms,3))
      allocate(fnnxc(natoms,3))
      allocate(fnnyc(natoms,3))


c  finding the first nearest neighbours
    do mind = 1,natoms
     
    fnnatm = 0
     do i = 1,3
       
        do j = 1,natoms

        xcatm = x(2,mind)
        ycatm = y(2,mind)
        
          xdist = x(i,j) - xcatm
          ydist = y(i,j) - ycatm
          distmag = xdist*xdist + ydist*ydist
          distmag = sqrt(distmag)
  
          if ((distmag < 1.42+eps ) .and. (distmag > 1.42-eps)) then
        fnnatm = fnnatm + 1
        fnn(mind,fnnatm) = j
        fnnxc(mind,fnnatm) = x(i,j)
        fnnyc(mind,fnnatm) = y(i,j)
        endif
  
        enddo
     enddo

    enddo
    
      
c generating points in the k space     
      
      lowlimit = -1.0*pi/a
      uplimit = 1.0*pi/a
      incremnt = 0.01
      
      vrange = uplimit-lowlimit
      numpoints = int(vrange/incremnt)

      print*, numpoints
      allocate (kx(numpoints))
      allocate (ky(numpoints))

      do i = 1,numpoints
      kx(i) = lowlimit + dble((i-1)*incremnt)
      ky(i) = lowlimit + dble((i-1)*incremnt)
      enddo
      
      
      nkx = size(kx);
      nky = size(ky);
      nkz = 1;


  
! generating the tight binding hamiltonian for a-gnr
      do mind = 1, size(kx)
    
      do nind = 1,size(ky)
      
        kxc = kx(mind);
        kyc = ky(nind);
        
c        print*, kxc,kyc
        
      do xx = 1,natoms
      do yy = 1,natoms
                
      H(xx,yy) = 0.0
                
      enddo
      enddo
        
         
      do i = 1,natoms
      do fnnj = 1,3
                   
      fnnatm = fnn(i,fnnj)
                
      if (fnnatm .ne. 0) then
                              
      numfnn = fnnyc(i,fnnj)- y(2,i)
      denomfnn = fnnxc(i,fnnj)- x(2,i)
                
               
      if (denomfnn > 0 .and. numfnn > 0) then
          
      H(i,fnnatm) = g1*exp( dcmplx(0,kxc*a1) + dcmplx(0,kyc*a2))

      elseif (denomfnn < 0 .and. numfnn < 0) then
      H(i,fnnatm) = g1*exp( dcmplx(0,-kxc*a1) + dcmplx(0,-kyc*a2))

      elseif (denomfnn < 0 .and. numfnn > 0) then
      H(i,fnnatm) = g1*exp(dcmplx(0,-kxc*a1) + dcmplx(0,kyc*a2))

      elseif (denomfnn > 0 .and. numfnn < 0) then
      H(i,fnnatm) = g1*exp(dcmplx(0,kxc*a1) + dcmplx(0,-kyc*a2)) 

      elseif (denomfnn > 0 .and. numfnn < 0 +eps .and. 
     +numfnn > 0-eps .and.denomfnn <1.42+eps .and. 
     +denomfnn > 1.42-eps) then
      H(i,fnnatm) = g1*exp(dcmplx(0,kxc*a3));

      elseif (denomfnn < 0 .and. numfnn < 0 +eps .and. 
     +numfnn > 0-eps .and. denomfnn < -1.42+eps .and. denomfnn 
     +> -1.42-eps) then
      H(i,fnnatm) = g1*exp(dcmplx(0,-kxc*a3));

      
      endif
                
      endif
      enddo
            
      
      
      H(i,i) = 0
      S(i,i) = 1
      enddo
                  
      call ZGEEV('N', 'N', natoms, H, natoms, b, DUMMY, 1, DUMMY,
     +1, WORK, 2*natoms, WORK, ok)


      write (252,10) (real(b(j)), j = 1, natoms), kxc,kyc
   10 format (16F10.3)
     
      enddo
      enddo


      close(252)
      close(253)
      
            
      
      end program
