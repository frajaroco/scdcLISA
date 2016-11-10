!
!     Francisco J. Rodriguez-Cortes, November 2016 
!
!     This function provides an edge-corrected kernel based estimator 
!     of the spatial second-order product density LISA functions given
!     in Cressie and Collins (2001).
!

       subroutine corelisa(i,d,n,s,ns,ks,hs,wrs,edge,lisa)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ks,ns,edge
       double precision kerns,s,hs,hij,wij,wrs,d,lisa
       dimension d(n,n),s(ns),wrs(n,n),lisa(ns),edge(2),ks(3)

       do iu=1,ns
        do j=1,n
          if (j.ne.i) then
           hij=d(i,j)
           if (ks(1).eq.1) then
               kerns=boxkernel((s(iu)-hij)/hs,hs)
                else if (ks(2).eq.1) then
                 kerns=ekernel((s(iu)-hij)/hs,hs)
                  else if (ks(3).eq.1) then
                   kerns=qkernel((s(iu)-hij)/hs,hs)
              end if
            if (kerns.ne.0d0) then
!     none   
              if (edge(1).eq.1) then
                 wij=kerns/s(iu)
                 lisa(iu)=lisa(iu)+wij 
              end if                         
!    isotropic
              if (edge(2).eq.1) then                  
                 wij=(kerns*wrs(i,j))/s(iu)
                 lisa(iu)=lisa(iu)+wij
              end if
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
