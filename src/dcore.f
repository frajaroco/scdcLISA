C
C     Francisco J. Rodriguez-Cortes, May 2017 
C
C     This function provides the coefficiebts and kernel basis of the 
C     spatial pair correlation LISA functions estimator.
C

       subroutine dcore(i,d,n,s,ns,ks,hs,lamd,wrs,edge,c,phi)

       implicit real*8(a-h,o-z)

       integer i,j,iu,n,ks,ns,edge,p
       double precision kerns,s,hs,hij,c,wrs,d,phi,lamd
       dimension d(n,n),s(ns),wrs(n,n)
       dimension c(n-1),phi((n-1),ns),edge(2),ks(3),lamd(n)
       
       do iu=1,ns
        p=0
        do j=1,n
          if (j.ne.i) then
           hij=d(i,j)
           p=p+1
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
                 c(p)=((n-1))/(lamd(i)*lamd(j))
                 phi(p,iu)=kerns/s(iu)
              end if                         
C    isotropic
              if (edge(2).eq.1) then  
                 c(p)=((n-1)*wrs(i,j))/(lamd(i)*lamd(j))
                 phi(p,iu)=kerns/s(iu)                
              end if
            end if
          end if
        end do
       end do

        return

        end
