!
!     Francisco J. Rodriguez-Cortes, November 2016 
!
!     This function provides an edge-corrected kernel based estimator 
!     of the spatial second-order product density function given in 
!     Cressie and Collins (2001).
!

       subroutine coresopdf(x,y,n,s,ns,ks,hs,area,wrs,edge,corepd)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ns,ks,edge
       double precision x,y,corepd,wij,s,hs,kerns
       double precision hij,two,xi,yi,pi,area
       dimension x(n),y(n),s(ns),wrs(n,n),edge(2),ks(3),corepd(ns)

          two=2d0
          pi=3.1415927d0

       do iu=1,ns
        do i=1,n
         xi=x(i)
         yi=y(i)      
          do j=1,n
           if (j.ne.i) then
            hij=sqrt((xi-x(j))**two+(yi-y(j))**two)
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
                 wij=kerns/(area*two*pi*s(iu))
                 corepd(iu)=corepd(iu)+wij 
              end if                         
!    isotropic
              if (edge(2).eq.1) then                  
                 wij=(kerns*wrs(i,j))/(area*two*pi*s(iu))
                 corepd(iu)=corepd(iu)+wij
              end if            
             end if
           end if
          end do
        end do
       end do
       
        return

        end
