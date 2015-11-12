!
!     Francisco J. Rodriguez-Cortes, November 2015
!
!     This function provides an edge-corrected kernel
!     estimator of the secon-order product density LISA function.
!

       subroutine corelisa(x,y,n,s,ns,ks,hs,xi,yi,area,i,wrs,lisa)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ks,ns
       real*8 kernes,x,y,s,hs,xj,yj,xi,yi,hij,wij,wrs,area,lisa
       dimension x(n),y(n),s(ns),wrs(n,n),lisa(ns)

       two=2d0
       pi=3.141592654d0

       do iu=1,ns
        lisa(iu)=0d0
        do j=1,n
         xj=x(j)
         yj=y(j)
          if (j.ne.i) then
           hij=sqrt((xi-xj)**two+(yi-yj)**two)
            if (ks.eq.1) then
             kerns=boxkernel((s(iu)-hij)/hs,hs)
              else if (ks.eq.2) then
               kerns=ekernel((s(iu)-hij)/hs,hs)
                else if (ks.eq.3) then
                 kerns=qkernel((s(iu)-hij)/hs,hs)
            end if
             if (kerns.ne.0d0) then
              wij=((n-1)*kerns*wrs(i,j))/(area*two*pi*s(iu))
              lisa(iu)=lisa(iu)+wij
            end if
          end if
        end do
       end do

        return

        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     functions called by :
!     -----------------------------------------
!
!     * boxkernel, ekernel, qkernel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------
!
!     boxkernel
!
!--------------------------------------------------------------------

       function boxkernel(x,h)

       implicit real*8 (a-h,o-z)

       double precision x, h

       if (abs(x).le.1d0) then
           boxkernel=1d0/2d0
       else
           boxkernel=0d0
       end if
       boxkernel=boxkernel/h

       return
       end

!--------------------------------------------------------------------
!
!     Epanechnikov kernel
!
!--------------------------------------------------------------------

       function ekernel(x,h)

       implicit real*8 (a-h,o-z)

       double precision x

       if (abs(x).le.1d0) then
           ekernel=(3d0/4d0)*(1-x**2)
       else
           ekernel=0d0
       end if
       ekernel=ekernel/h

       return
       end

!--------------------------------------------------------------------
!
!     quartic (biweight) kernel
!
!--------------------------------------------------------------------

       function qkernel(x,h)

       implicit real*8 (a-h,o-z)

       double precision x, h

       if (abs(x).le.1d0) then
           qkernel=(15d0/16d0)*(1-x**2)**2
       else
           qkernel=0d0
       end if
       qkernel=qkernel/h

       return
       end
