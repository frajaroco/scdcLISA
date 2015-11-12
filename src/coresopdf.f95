!
!     Francisco J. Rodriguez-Cortes, November 2015 
!
!     This function provides an edge-corrected kernel 
!     estimator of the spatial secon-order product density function.
!

       subroutine coresopdf(x,y,n,s,ns,ks,hs,area,wrs,corepd)   	   

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ns,ks
       double precision x,y,corepd,wij,s,hs,kerns
       double precision hij,two,xi,yi,pi,area
       dimension x(n),y(n),s(ns),wrs(n,n),corepd(ns)

          two=2d0
          pi=3.1415927d0
		 
       do iu=1,ns
        do i=1,n
         xi=x(i)
         yi=y(i)      
          do j=1,n
           if (j.ne.i) then
            hij=sqrt((xi-x(j))**two+(yi-y(j))**two)
            if (ks.eq.1) then
             kerns=boxkernel((s(iu)-hij)/hs,hs)
              else if (ks.eq.2) then
               kerns=ekernel((s(iu)-hij)/hs,hs)
                else if (ks.eq.3) then
                 kerns=qkernel((s(iu)-hij)/hs,hs)
            end if   
             if (kerns.ne.0d0) then
              wij=(kerns*wrs(i,j))/(area*two*pi*s(iu))
              corepd(iu)=corepd(iu)+wij
			 end if
           end if
          end do
        end do
       end do
       
        return

        end
        
