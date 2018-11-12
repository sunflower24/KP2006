subroutine guess2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

use params
implicit none

! fills the autarky functions with permanent income consumption

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do rc=1,pm

			cmfun(tc,tyc,sc,ic,rc)=((Rm(tc)-1.0)/Rm(tc))*( gridm(tyc,rc)+expinc(tc,tyc,sc,ic) )
			if ( Rm(tc) < 1.00001) then
				cmfun(tc,tyc,sc,ic,rc)=((Rm(tc)-0.9)/Rm(tc))*( gridm(tyc,rc)+expinc(tc,tyc,sc,ic) )
			endif
			sfun(tc,tyc,sc,ic,rc)=Rm(tc)* ( w(tc)*yat(tc,tyc,sc,ic)+gridm(tyc,rc)-cmfun(tc,tyc,sc,ic,rc) )
			if ( sfun(tc,tyc,sc,ic,rc) < 0.0 ) then
				sfun(tc,tyc,sc,ic,rc)=0.0
				cmfun(tc,tyc,sc,ic,rc)=w(tc)*yat(tc,tyc,sc,ic)+gridm(tyc,rc)
			end if 
			vmfun(tc,tyc,sc,ic,rc)=U(cmfun(tc,tyc,sc,ic,rc))
			vpmfun(tc,tyc,sc,ic,rc)=(1.0-bta)*MU(cmfun(tc,tyc,sc,ic,rc))
			end do
		end do
	end do
end do



end subroutine guess2