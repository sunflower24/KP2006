subroutine autval

 use params

 implicit none
 real(prec),dimension(ns*ni,ns*ni)::eye,mat,matinv
 real(prec),dimension(ns*ni)::put,helper
 integer i

do tyc=1,ntypes
	eye=0.0
	do sc=1,ns
		do ic=1,ni
			put((sc-1)*ni+ic)=(1.0-bta)*U( w(tc)*yat(tc,tyc,sc,ic) )
			eye((sc-1)*ni+ic,(sc-1)*ni+ic)=1.0
			do scc=1,ns
				do icc=1,ni
					mat((sc-1)*ni+ic,(scc-1)*ni+icc)=eye((sc-1)*ni+ic,(scc-1)*ni+icc)-bta*probs(tc,sc,scc)*pi(tc,icc)
				end do
			end do 
		end do
	end do
	call dlinrg(ns*ni,mat,ns*ni,matinv,ns*ni)
	helper=matmul(matinv,put)
	do sc=1,ns
		vmaut(tc,tyc,sc,1:ni)=helper((sc-1)*ni+1:sc*ni)
	end do

	
	eye=0.0
	do sc=1,ns
		do ic=1,ni
			put((sc-1)*ni+ic)=w(tc)*yat(tc,tyc,sc,ic)
			eye((sc-1)*ni+ic,(sc-1)*ni+ic)=1.0
			do scc=1,ns
				do icc=1,ni
					mat((sc-1)*ni+ic,(scc-1)*ni+icc)=eye((sc-1)*ni+ic,(scc-1)*ni+icc)-(1.0/R(tc))*probs(tc,sc,scc)*pi(tc,icc)
				end do
			end do 
		end do
	end do
	call dlinrg(ns*ni,mat,ns*ni,matinv,ns*ni)
	helper=matmul(matinv,put)
	do sc=1,ns
		expinc(tc,tyc,sc,1:ni)=helper((sc-1)*ni+1:sc*ni)
	end do
end do
 
if (choice2>1) then
	call autvalsav
end if


end subroutine autval

 