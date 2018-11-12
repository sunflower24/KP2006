subroutine guess                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

use params
implicit none

! fills the function with the complete markets solution	

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa

				cfun(tc,tyc,sc,ic,ac)=(1.0- ( (R(tc)*bta)**(1.0/sigma) )/R(tc) )*( grida(tc,tyc,sc,ic,ac)+expinc(tc,tyc,sc,ic) )
				if (abs(sigma-1.0) >0.0001) then
					vfun(tc,tyc,sc,ic,ac)=( (1.0-bta)/(1.0-bta*(bta*R(tc))**((1.0-sigma)/sigma) ) )*U(cfun(tc,tyc,sc,ic,ac))
				else
					vfun(tc,tyc,sc,ic,ac)=U(cfun(tc,tyc,sc,ic,ac))+(bta/(1.0-bta))*log(bta*R(tc))
				endif
				lafun(tc,tyc,sc,ic,ac)=(1.0-bta)*MU(cfun(tc,tyc,sc,ic,ac))
				do scc=1,ns
					do icc=1,ni
						afun(tc,tyc,sc,ic,ac,scc,icc)= ( (bta*R(tc))**(1.0/sigma) )*( grida(tc,tyc,sc,ic,ac)+expinc(tc,tyc,sc,ic))-expinc(tc,tyc,scc,icc) 
					end do
				end do
			end do
		end do
	end do
end do


end subroutine guess