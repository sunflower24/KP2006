subroutine guesstr                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

use params
implicit none

! fills the function with optimal policy from final SS

if (tc>2) then
  do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa

			cfun(tc,tyc,sc,ic,ac)=cfun(nt,tyc,sc,ic,ac)
			vfun(tc,tyc,sc,ic,ac)=vfun(nt,tyc,sc,ic,ac)
			lafun(tc,tyc,sc,ic,ac)=lafun(nt,tyc,sc,ic,ac)
	
				do scc=1,ns
					do icc=1,ni
						afun(tc,tyc,sc,ic,ac,scc,icc)=afun(nt,tyc,sc,ic,ac,scc,icc)
					end do
				end do
			end do
		end do
	end do
  end do
end if

if (tc==2) then
  do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa

			cfun(tc,tyc,sc,ic,ac)=cfun(1,tyc,sc,ic,ac)
			vfun(tc,tyc,sc,ic,ac)=vfun(1,tyc,sc,ic,ac)
			lafun(tc,tyc,sc,ic,ac)=lafun(1,tyc,sc,ic,ac)
	
				do scc=1,ns
					do icc=1,ni
						afun(tc,tyc,sc,ic,ac,scc,icc)=afun(1,tyc,sc,ic,ac,scc,icc)
					end do
				end do
			end do
		end do
	end do
  end do
end if


end subroutine guesstr