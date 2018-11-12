subroutine saving

use params
implicit none
integer::cc
real(prec)::giniy,ginic
real(prec),dimension(ni*ns*pa*ntypes)::lorincx,lorincy,lorconsx,lorconsy
real(prec),dimension(nt,ni*ns*pa*ntypes)::lorenzyx,lorenzyy,lorenzcx,lorenzcy

open(unit=98,file='modeldist.txt')

! Restore initial and final distribution

tc=nt
omega(nt,1:ntypes,1:ns,1:ni,1:pa)=finomega
call gini2(omega(nt,1:ntypes,1:ns,1:ni,1:pa),lorincx,lorincy,lorconsx,lorconsy,giniy,ginic)

! Save Distributions

write(98,*) 'Distributions of Assets'

do tc=1,nt
  do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				write(98,'(2f9.6)') grida(tc,tyc,sc,ic,ac), omega(tc,tyc,sc,ic,ac)
			end do
		end do
	end do
  end do
end do 

write(98,*) 'Distributions of Income'

do tc=1,nt
	do tyc=1,ntypes
		do sc=1,ns
			do ic=1,ni
				write(98,'(2f9.6)') w(tc)*yat(tc,tyc,sc,ic), diste(tc,sc)*pi(tc,ic)
			end do
		end do
	end do
end do 

write(98,*) 'Distributions of Consumption'

do tc=1,nt
	do tyc=1,ntypes
		do sc=1,ns
			do ic=1,ni
				do cc=1,pc
					write(98,'(2f9.6)') gridc(cc), cdist(tc,tyc,cc,sc,ic)
				end do
			end do
		end do
	end do
end do 


! Compute and save total intermediation

write(98,*) 'Total Consumer Credit'
totint=0.0
do tc=1,nt
write(98,'(i4,f9.6)') tc,totint(tc)

end do 

! Save Gini Coefficients and other Dispersion Measures of Income and Consumption

write(98,*) 'Std dev of log income and log consumption'

do tc=1,nt
	write(98,'(2f9.6)') stdly(tc),stdlc(tc)
end do 


write(98,*) 'Gini Coefficients and Lorenz Curves'

do tc=1,nt
	call gini2(omega(tc,1:ntypes,1:ns,1:ni,1:pa),lorincx,lorincy,lorconsx,lorconsy,giniy,ginic)
	giniyvec(tc)=giniy
	ginicvec(tc)=ginic
	lorenzyx(tc,1:pa*ni*ns*ntypes)=lorincx
	lorenzyy(tc,1:pa*ni*ns*ntypes)=lorincy
	lorenzcx(tc,1:pa*ni*ns*ntypes)=lorconsx
	lorenzcy(tc,1:pa*ni*ns*ntypes)=lorconsy
end do
write(98,*) 'Gini Coefficients for Income and Consumption'
do tc=1,nt
	write(98,'(2f9.6)') ginicvec(tc), giniyvec(tc)
end do
write(98,*) 'Initial and final Lorenz Curves for Income and Consumption'
tc=1
	do ac=1,pa*ni*ns*ntypes
		write(98,'(4f9.6)') lorenzyx(tc,ac),lorenzyy(tc,ac),lorenzcx(tc,ac),lorenzcy(tc,ac)
	end do

tc=nt
	do ac=1,pa*ns*ntypes
		write(98,'(4f9.6)') lorenzyx(tc,ac),lorenzyy(tc,ac),lorenzcx(tc,ac),lorenzcy(tc,ac)
	end do



! Save Equilibrium Interest Rates

write(98,*) 'Equilibrium Interest Rates'

do tc=1,nt

write(98,'(1f9.6)') R(tc)

end do

close(98)

end subroutine saving