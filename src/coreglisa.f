C
C     Francisco J. Rodriguez-Cortes, March 2017 
C
C     This function provides an edge-corrected kernel based estimator 
C     of the spatial pair correlation LISA functions.
C

       subroutine coreglisa(i,d,n,s,ns,ks,hs,lamd,wrs,edge,glisa)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ks,ns,edge
       double precision kerns,s,hs,hij,wij,wrs,d,glisa,lamd
       dimension d(n,n),s(ns),wrs(n,n),glisa(ns),edge(2),ks(3),lamd(n)

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
C     none   
              if (edge(1).eq.1) then
                 wij=kerns/(lamd(i)*lamd(j)*s(iu))
                 glisa(iu)=glisa(iu)+wij 
              end if                         
C    isotropic
              if (edge(2).eq.1) then                  
                 wij=(kerns*wrs(i,j))/(lamd(i)*lamd(j)*s(iu))
                 glisa(iu)=glisa(iu)+wij
              end if
            end if
          end if
        end do
       end do

        return

        end
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C
C     functions called by :
C     -----------------------------------------
C
C     * boxkernel, ekernel, qkernel
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C--------------------------------------------------------------------
C
C     boxkernel
C
C--------------------------------------------------------------------

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

C--------------------------------------------------------------------
C
C     Epanechnikov kernel
C
C--------------------------------------------------------------------

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

C--------------------------------------------------------------------
C
C     quartic (biweight) kernel
C
C--------------------------------------------------------------------

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