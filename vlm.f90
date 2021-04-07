      program vlm
!     
!     VORTEX LATTICE METHOD     
!     
      implicit none 
      integer      nx,ny
      parameter   (nx=101)
      parameter   (ny=101)
      integer      Np,Np_w
      parameter   (Np=10201)
      parameter   (Np_w=(50)*10201)
      real*8       hlx, hly     
      real*8       xn(nx+1,ny+1,3)
      real*8       xn_w(nx+1,(50)*ny+1,3)
      real*8       xp(3,Np,5)
      real*8       xp_w(3,Np_w,5)
      real*8       xm(Np,4,3)
      real*8       xc(3,Np)
      real*8       r(3,Np,Np,4)
      real*8       xpl(Np,4)
      real*8       Uinf(3)
      real*8       vers_incl(3,Np)
      real*8       versp(3,Np) 
      real*8       alpha
       
      call input(nx,ny,hlx,hly,alpha)
      call node_geom(nx,ny,hlx,hly)
      call node_geom_wake(nx,ny,hlx,hly)
      call out_geom(nx,ny,xn)
      call vert_panel(nx,ny,xn)
      call vert_panel_wake(nx,ny,xn_w)
      call center_geom(xp)
      call vers_panel(xp)
      call incl_vers(alpha,versp)      
      call side_panel(xp)
      call side_panel_wake(xp_w)
      call vect_r(xp,xc)
      call term_noti(Uinf,vers_incl)
      call influence_matrix(xpl,r,xc,xm)
      
      end program
      
 
      subroutine input(nx,ny,hlx,hly,alpha)
!     
!     nx       = discretizzazione lungo x
!     ny       = discretizzazione lungo y 
!     hlx      = apertura alare
!     hly      = corda
!     alpha    = angolo d'attacco
!     Uinf_sca = intensità della corrente indisturbata
      
      implicit none
      
      integer nx,ny
      real*8  hlx, hly
      real*8  alpha
      real*8  Uinf(3)
      real*8  Uinf_sca
      
      open(unit=1,file='vlm_input.dat')
         read(1,*) nx	  
         read(1,*) ny	  
         read(1,*) hlx	  
         read(1,*) hly
         read(1,*) alpha
         read(1,*) Uinf_sca
      close(1)
      
      Uinf(1) = 0.d0
      Uinf(2) = Uinf_sca
      Uinf(3) = 0.d0
      
      return
      end subroutine	  

              
      subroutine node_geom(nx,ny,hlx,hly)
!     
!     xn(i,j,k)= coordinate dei nodi dell'ala
!     k=1,2,3
!     i pannelli dell'ala si numerano da sinistra verso destra
!     e dal basso verso l'alto
!     
      implicit none
      real*8  Dx, Dy
      real*8  hlx, hly
      integer nx, ny
      integer i,j
      real*8  xx,yy
      real*8  xn(nx+1,ny+1,3) 
      
      Dx=hlx/nx
      Dy=hly/ny
      
!     Dx = lunghezza del singolo  pannello lungo x
!     Dy = lunghezza del singolo  pannello lungo y
!     
      do i=1, nx+1
         xx=float(i-1)*Dx
         do j=1, ny+1
            yy=float(j-1)*Dy
            xn(i,j,1)=xx
            xn(i,j,2)=yy
            xn(i,j,3)=0.d0
         end do    
      end do
      
      end subroutine

      
      subroutine node_geom_wake(nx,ny,hlx,hly)
!     
!     xn(i,j,k)= coordinate dei nodi della scia
!     k=1,2,3
!     i pannelli della scia si numerano da sinstra verso destra e 
!     dall'alto verso il basso 
!           
      implicit none
      real*8  Dx, Dy
      real*8  hlx, hly
      integer nx, ny
      integer i,j
      real*8  xx,yy
      real*8  xn_w(nx+1,(50)*ny+1,3)       
      
      Dx=hlx/nx
      Dy=hly/ny
!      
!     Dx = lunghezza del singolo  pannello lungo x
!     Dy = lunghezza del singolo  pannello lungo y
!     
      do i=1, nx+1
         xx=float(i-1)*Dx
         do j=1,(50*ny)+1
            yy=-float(j-1)*Dy
            xn_w(i,j,1)=xx
            xn_w(i,j,2)=yy
            xn_w(i,j,3)=0.d0
         end do    
      end do
      
      end subroutine
 
     
      subroutine out_geom(nx,ny,xn)
!     
!     
!     
      
      implicit none
      integer    i, j
      integer    nx, ny
      real*8     xx, yy, zz
      real*8     xn(nx+1,ny+1,3)    
      
      open(unit=1, file= 'wing_geometry.dat')
         write(1,*) nx,ny
         do i=1,nx+1
            do j=1,ny+1
               xx=xn(i,j,1) 
               yy=xn(i,j,2) 
               zz=xn(i,j,3) 
               write(1,*) xx, yy, zz
            end do
         end do
      close(1)
      
      end subroutine
      

      subroutine vert_panel(nx,ny,xn)
!      
!     xp(k,p,:) individua i 5 vertici intorno al p-esimo pannello
!         
      implicit none
      integer    p,k
      integer    i,j
      integer    Np
      integer    nx,ny
      parameter   (Np=10201)
      real*8     xp(3,Np,5)
      real*8     xn(nx+1,ny+1,3)
       
      do p=1,Np
         do k=1,3
            if (p<nx) then
               i=p
               j=0 
               else if ((mod(p,nx))==0) then
                  i=nx
                  j=p/nx
                  else
                     i=mod(p,nx)
                     j=(p/nx)+1
                     xp(k,p,1)=xn(i,j,k)
                     xp(k,p,2)=xn(i+1,j,k)
                     xp(k,p,3)=xn(i+1,j+1,k)
                     xp(k,p,4)=xn(i,j+1,k)
                     xp(k,p,5)=xp(k,p,1)             
            end if 
         end do 
      end do       
      
      end subroutine
      
      
      subroutine vert_panel_wake(nx,ny,xn_w)
!     
!     xp_w(k,p,:) individua i 5 vertici intorno al p-esimo pannello della scia
!         
      implicit none
      integer    p,k
      integer    i,j
      integer    Np
      integer    Np_w
      integer    nx,ny
      parameter  (Np_w=(50)*10201)
      parameter  (Np=10201)
      real*8     xp_w(3,Np_w,5)
      real*8     xn_w(nx+1,(50)*ny+1,3)      
      
      do p=1,Np_w
         do k=1,3
            if (p<nx) then
               i=p
               j=0 
               else if ((mod(p,nx))==0) then
                  i=nx
                  j=p/nx
                  else
                     i=mod(p,nx)
                     j=(p/nx)+1
                     xp_w(k,p,1)=xn_w(i,j,k)
                     xp_w(k,p,2)=xn_w(i+1,j,k)
                     xp_w(k,p,3)=xn_w(i+1,j+1,k)
                     xp_w(k,p,4)=xn_w(i,j+1,k)
                     xp_w(k,p,5)=xp_w(k,p,1)             
            end if 
         end do 
      end do
      
      end subroutine       
      
      
      subroutine center_geom(xp)
!     
!     xc(:,p) contiene le coordinate dei centri dei pannelli
!           
      implicit none
      integer    p
      integer    Np
      parameter  (Np=10201)
      real*8     xp(3,Np,5)
      real*8     xc(3,Np)
      
      do p=1,Np
         xc(1,p)=(xp(1,p,1)+xp(1,p,2))/2
         xc(2,p)=(xp(2,p,1)+xp(2,p,4))/2
         xc(3,p)=0.d0
      end do
      
      end subroutine
      
      
      subroutine center_geom_wake(xp_w)
!     
!     xc(:,p) contiene le coordinate dei centri dei pannelli della scia
!     
      implicit none
      integer    p
      integer    Np_w
      parameter  (Np_w=(50)*10201)
      real*8     xp_w(3,Np_w,5)
      real*8     xc_w(3,Np_w)
      
      do p=1,Np_w
         xc_w(1,p)=(xp_w(1,p,1)+xp_w(1,p,2))/2
         xc_w(2,p)=(xp_w(2,p,1)+xp_w(2,p,4))/2
         xc_w(3,p)=0.d0
      end do
      
      end subroutine
      
      
      subroutine vers_panel(xp)
!     
!     vers(p) = versore normale al p-esimo pannello
!           
      implicit none
      real*8     v1(3), v2(3) 
      integer    p
      integer    Np
      parameter  (Np=10201)
      real*8     xp(3,Np,5)
      real*8     norm_crossp
      real*8     versp(3,Np) 
      real*8     module      
      real*8     crossp(3)
       
      do p=1,Np
         v1(1)=xp(1,p,3)-xp(1,p,1)
         v1(2)=xp(2,p,3)-xp(2,p,1)
         v1(3)=xp(3,p,3)-xp(3,p,1)
         v2(1)=xp(1,p,2)-xp(1,p,4)
         v2(2)=xp(2,p,2)-xp(2,p,4)
         v2(3)=xp(3,p,2)-xp(3,p,4)
          
!        calcolo il prodotto vettoriale e la norma
      
         crossp(1)=v1(2)*v2(3)-v1(3)*v2(2)
         crossp(2)=v1(3)*v2(1)-v1(1)*v2(3)
         crossp(3)=v1(1)*v2(2)-v1(2)*v2(1)
 
         norm_crossp=module(crossp)
       
         versp(1,p)=crossp(1)/norm_crossp
         versp(2,p)=crossp(2)/norm_crossp
         versp(3,p)=crossp(3)/norm_crossp   
      
      end do
      
      end subroutine   
      
      
      subroutine incl_vers(alpha,versp)
!     
!     vers_incl(:,p) =  vettore normale al p-esimo pannello inclinato di alpha gradi
!     
      implicit none
      integer     p, i, j
      integer     Np
      parameter  (Np=10201)
      real*8      R1(3,3)
      real*8      vv(3)
      real*8      versp(3,Np) 
      real*8      vers_incl(3,Np)
      real*8      alpha
      real*8      sa
      real*8      ca
      real*8      pi
      parameter   (pi=3.141592654)
      
      alpha=alpha*pi/180
      ca=cos(alpha)
      sa=sin(alpha)
      
!     R1 = matrice di rotazione attorno al bordo d'attacco                  
      
      R1(1,1)=1      
      R1(1,2)=0      
      R1(1,3)=0      
      R1(2,1)=0      
      R1(2,2)=ca      
      R1(2,3)=sa     
      R1(3,1)=0      
      R1(3,2)=-sa      
      R1(3,3)=ca      
      
      do p=1,Np
         do i=1,3
            do j=1,3
               vv(j)=R1(i,j)*versp(j,p)
            end do
            vers_incl(i,p)=vv(1)+vv(2)+vv(3)
         end do
      end do
       
      end subroutine
      
      
      subroutine side_panel(xp)
!     
!     xpl(p,l) = lunghezza del l-esimo lato del p-esimo pannello 
!     
      implicit none
      integer     p,l,k
      integer     Np
      parameter  (Np=10201)
      real*8      xd1(3)      
      real*8      xd2(3)      
      real*8      xp(3,Np,5)
      real*8      xpl(Np,4)
      real*8      distanza
      
      do p=1,Np
         do l=1,4
            do k=1,3
               xd1(k)=xp(k,p,l+1)
               xd2(k)=xp(k,p,l)
            end do
            xpl(p,l)=distanza(xd1,xd2)          
         end do
      end do     
                 
      end subroutine   
      
                       
      subroutine side_panel_wake(xp_w)
!                     
!     xpl_w(p,l) = lunghezza del l-esimo lato del p-esimo pannello della scia 
!                     
      implicit none    
      integer     p,l,k
      integer     Np_w 
      parameter  (Np_w=(50)*10201)
      real*8      xd1_w(3)      
      real*8      xd2_w(3)      
      real*8      xp_w(3,Np_w,5)
      real*8      xpl_w(Np_w,4)
      real*8      distanza
                      
      do p=1,Np_w     
         do l=1,4     
            do k=1,3   
               xd1_w(k)=xp_w(k,p,l+1)
               xd2_w(k)=xp_w(k,p,l)
            end do    
            xpl_w(p,l)=distanza(xd1_w,xd2_w)          
         end do      
      end do         
      
      end subroutine
      
      
      subroutine term_noti(Uinf,vers_incl)
!     
!     T = vettore dei termini noti 
!     
      implicit none
      integer     Np,p
      parameter  (Np=10201)
      real*8      T(Np)
      real*8      Uinf(3)
      real*8      vers_incl(3,Np)
      real*8      vn(3)
      real*8      prod
      
      do p=1,Np
          vn(1)=vers_incl(1,p)
          vn(2)=vers_incl(2,p)
          vn(3)=vers_incl(3,p)
          T(p)=-prod(Uinf,vn) 
      end do      
      
      end subroutine    
     
      subroutine punto_medio(xp)
!     
!     xm(p,l,:) = punto medio del' l-esimo lado del p-esimo pannello 
!     
      implicit none
      integer     p,l,k
      integer     Np
      parameter  (Np=10201)
      real*8      xm(Np,4,3)
      real*8      xp(3,Np,5)
      
      do p=1,Np
         do l=1,4
            do k=1,3
               xm(p,l,k)=(xp(k,p,l)+xp(k,p,l+1))/2
            end do
         end do
      end do      
      
      end subroutine
      
      
      subroutine punto_medio_wake(xp_w)
!     
!     xm_w(p,l,:) = punto medio del' l-esimo lado del p-esimo pannello della scia
!     
      implicit none
      integer     p,l,k
      integer     Np_w
      parameter  (Np_w=(50)*10201)
      real*8      xm_w(Np_w,4,3)
      real*8      xp_w(3,Np_w,5)
      
      do p=1,Np_w
         do l=1,4
            do k=1,3
               xm_w(p,l,k)=(xp_w(k,p,l)+xp_w(k,p,l+1))/2
            end do
         end do
      end do      
      
      end subroutine


      subroutine vect_r(xm,xc)
     
     r(:,p,l) è il vettore che individua il punto di controllo del p-esimo pannello dal centro del filamento 
     del j-esimo pannello  

      implicit none
      integer    p,l,k,j
      integer    Np
      parameter (Np=10201)
      real*8     xm(Np,4,3)
      real*8     xc(3,Np)
      real*8     r(3,Np,Np,4)
      
      do p=1,Np
         do j=1,Np
            do l=1,4
               do k=1,3
                  r(k,p,j,l)=xm(j,l,k)-xc(k,p)
               end do
            end do
         end do
      end do   
      
      end subroutine      
      
      
      subroutine influence_matrix(xpl,r,xc,xm)
     
!     AIC = matrice di influenza
!     AIC(p,l) contiene il valore della velocità indotta nell' l-esimo  punto di controllo
!            dal vortice unitario posto nel  p-esimo pannello 
!     dl e r sono ortogonali
      
      implicit none
      integer      p,l,k,j
      integer      Np
      parameter   (Np=10201)
      real*8       pi
      parameter   (pi=3.141592654)
      real*8       rr
      real*8       dlk(4)
      real*8       xm(Np,4,3)
      real*8       xc(3,Np)
      real*8       xpl(Np,4)
      real*8       rk(3)
      real*8       r(3,Np,Np,4)
      real*8       module
      real*8       Vind(4)
      real*8       VAPJ
      real*8       steta
      real*8       AIC(Np,Np)
      real*8       dx,dy       

      do p=1,Np
         do j=1,Np
             do l=1,4 
                dlk(l)=xpl(p,l)  
                   do k=1,3
                      rk(k)=xm(j,l,k)-xc(k,p)
                   end do               
                   rr=module(rk)
                   dy=abs(xm(j,l,2)-xc(2,p))
                   dx=abs(xm(j,l,1)-xc(1,p))
                   if ((mod(l,2))==0) then
                      steta=dx/rr
                      else
                         steta=dy/rr 
                   end if
                   Vind(l)=(dlk(l)/(4*pi*(rr**2)))*steta
              end do
              VAPJ=sum(Vind)
         end do 
         AIC(p,j)=VAPJ         
      end do
      
      end subroutine 
      
      
      real*8 function distanza(x,y)
!     
!     la funzione calcola la distanza tra due punti
!     
      implicit none
      real*8  d1,d2,d3
      real*8  x(3),y(3)
      real*8  z(3)
      real*8  module
      
      z=x-y
      
      distanza=module(z)
      
      return
      
      end function distanza
      
      
      real*8 function module(x)
!     
!     la funzione calcola il modulo di un vettore
!     
      implicit none 
      real*8  x(3)
      
      module=sqrt(x(1)**2+x(2)**2+x(3)**2)
      
      return
      
      end function module
      
      
      real*8 function prod(x,y)
!     
!     la funzione calcola il prodotto scalare tra due vettori
!     
      implicit none
      integer   i
      real*8    x(3), y(3)
                 
      do i=1,3
         prod=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
      end do
       
      return
      
      end function prod  
      
      


      
      
      
      
      
      


          






























