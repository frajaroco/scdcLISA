!
!     Francisco J. Rodriguez-Cortes, November 2016 
!
!     This function provides an edge-corrected kernel based estimator 
!     of the spatial second-order product density function given in 
!     Cressie and Collins (2001).
!

       subroutine coresopdf(d,n,s,ns,ks,hs,wrs,edge,corepd)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ns,ks,edge
       double precision corepd,hij,wij,s,kerns,d,hs
       dimension d(n,n),s(ns),wrs(n,n),edge(2),ks(3),corepd(ns)
        
        corepd=0d0

       do iu=1,ns
        do i=1,n 
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
                 corepd(iu)=corepd(iu)+wij 
              end if                         
!    isotropic
              if (edge(2).eq.1) then                  
                 wij=(kerns*wrs(i,j))/s(iu)
                 corepd(iu)=corepd(iu)+wij
              end if
             end if 
            end if           
          end do
        end do
       end do
       
        return

        end
