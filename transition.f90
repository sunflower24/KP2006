! This is a subroutine to compute the transition path

subroutine transition

use params
implicit none


integer::i,maxit=1000,upd
real(prec)::diff,adj=0.05
real(prec),dimension(nt)::Knew
real(prec), dimension(nt,ntypes,ns,ni,pa)::oldgrida,oldcfun,oldvfun,oldlafun
real(prec), dimension(nt,ntypes,ns,ni,pa,ns,ni)::oldafun
real(prec)::vals(2)
integer::inds(2)
real(prec),dimension(ns)::sol1,gue1

external borrowtr


print*,'Transition Started'

! Guess Initial Interest rates

help=amin
grida(2,:,:,:,:)=grida(1,:,:,:,:)
do tc=2,nt-1
	R(tc)=R(nt)
	amin(tc,:,:,:)=amin(nt,:,:,:)
	call dgrid
	call guesstr
end do

do i=1,maxit
print*,i

do tc=2,nt-1
	KA(tc)=LS(tc)*(A*alpha/(R(tc-1)-1.0+delta))**(1.0/(1.0-alpha))
	w(tc)=MPL(KA(tc),LS(tc))
end do


! compute optimal policy functions

	do ttc=2,nt-1
		tc=nt+1-ttc

		call autvaltr
 		call enforcetr

	end do

! compute cross-sectional distributions


	call disstr


! Check for Convergence

diff=maxval(abs(resd1))
print*,diff

if (diff<tol*0.01) then
	exit
end if


! Update guesses for R
Rnew(1)=R(1)
Rnew(nt)=R(nt)

do tc=2,nt
	Knew(tc)=assdemand(tc)/R(tc-1)
end do

do tc=3,nt
	Rnew(tc-1)=MPC(Knew(tc),LS(tc))+1.0
end do


open(unit=90,file='interest.txt')

do tc=1,nt
	write(90,'(2f9.6)') R(tc), Rnew(tc)
end do

close(90)

do tc=1,nt
	R(tc)=adj*Rnew(tc)+(1.0-adj)*R(tc)
end do

open(unit=20,file='GiniCoefficient.txt')
do tc=1,nt
	write(20,'(3f9.4)') giniyvec(tc), ginicvec(tc), totint(tc)
end do
close(20)

open(unit=22,file='Constraints.txt')
do tc=1,nt
	write(22,'(18f9.4)') amin(tc,:,:,:)
end do
close(22)

open(unit=24,file='Stddev.txt')
do tc=1,nt
	write(24,'(10f9.4)') varvec(tc,1:4), stdsvec(tc,1:4), varlc(tc), stdlc(tc)
end do
close(24)


! Update asset grid

do tc=2,nt-1

	amin(tc,:,:,:)=help(tc,:,:,:)

 	call dgrid

end do

! Redefine all functions to lie on the new grid


do tc=2,nt-1
	oldgrida(tc,:,:,:,:)=grida(tc,:,:,:,:)
	oldcfun(tc,:,:,:,:)=cfun(tc,:,:,:,:)
	oldafun(tc,:,:,:,:,:,:)=afun(tc,:,:,:,:,:,:)
	oldvfun(tc,:,:,:,:)=vfun(tc,:,:,:,:)
	oldlafun(tc,:,:,:,:)=lafun(tc,:,:,:,:)
end do

do tc=2,nt-1
  do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				call basefun(oldgrida(tc,tyc,sc,ic,:),pa,grida(tc,tyc,sc,ic,ac),vals,inds)
				cfun(tc,tyc,sc,ic,ac)=vals(1)*oldcfun(tc,tyc,sc,ic,inds(1))+vals(2)*oldcfun(tc,tyc,sc,ic,inds(2))
				vfun(tc,tyc,sc,ic,ac)=vals(1)*oldvfun(tc,tyc,sc,ic,inds(1))+vals(2)*oldvfun(tc,tyc,sc,ic,inds(2))
				lafun(tc,tyc,sc,ic,ac)=vals(1)*oldlafun(tc,tyc,sc,ic,inds(1))+vals(2)*oldlafun(tc,tyc,sc,ic,inds(2))
				do scc=1,ns
					do icc=1,ni
						afun(tc,tyc,sc,ic,ac,scc,icc)=vals(1)*oldafun(tc,tyc,sc,ic,inds(1),scc,icc)+vals(2)*oldafun(tc,tyc,sc,ic,inds(2),scc,icc)
					end do
				end do
			end do
		end do
	end do
  end do
end do

end do

end subroutine transition
