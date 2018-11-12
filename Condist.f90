subroutine condist		! Compute the distribution of consumption


use params
implicit none
integer inds(2),cc
real(prec),dimension(2)::vals
real(prec) mass,testc

	open(unit=16,file='cdist.txt')


cdist(tc,1:ntypes,1:pc,1:ns,1:ni)=0.0

do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			do ac=1,pa
				mass=omega(tc,tyc,sc,ic,ac)
				call basefun(gridc,pc,cfun(tc,tyc,sc,ic,ac),vals,inds)

				if ( (minval(vals)<0.0) .and. ( mass>0.00001 ) ) then
					print*,'Enlarge your consumption grid'
					print*,'Consumption is ',cfun(tc,tyc,sc,ic,ac)
				end if

				cdist(tc,tyc,inds(1),sc,ic)=cdist(tc,tyc,inds(1),sc,ic)+vals(1)*mass
				cdist(tc,tyc,inds(2),sc,ic)=cdist(tc,tyc,inds(2),sc,ic)+vals(2)*mass
			end do
		end do
	end do
end do

exput1(tc)=0.0

do tyc=1,ntypes
	do cc=1,pc
		do scc=1,ns
			do icc=1,ni
				exput1(tc)=exput1(tc)+U(gridc(cc))*cdist(tc,tyc,cc,scc,icc)
			end do
		end do
	end do
end do


rewind(16)
write(16,fmt=*) cdist,gridc  
!print*,' Cons. Distr. saved '
rewind(16)

! Computing standard deviation of log consumption

avglc(tc)=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			avglc(tc)=avglc(tc)+sum(log(gridc)*cdist(tc,tyc,1:pc,sc,ic))
		end do
	end do
end do

varlc(tc)=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			varlc(tc)=varlc(tc)+sum(((log(gridc)-avglc(tc))**2)*cdist(tc,tyc,1:pc,sc,ic))
		end do
	end do
end do
stdlc(tc)=sqrt(varlc(tc))

! By groups

tavglc(tc,:)=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			tavglc(tc,tyc)=tavglc(tc,tyc)+sum(log(gridc(1:pc))*(cdist(tc,tyc,1:pc,sc,ic)/pitypes(tyc)))
		end do
	end do
end do

tvarlc(tc,:)=0.0
do tyc=1,ntypes
	do sc=1,ns
		do ic=1,ni
			tvarlc(tc,tyc)=tvarlc(tc,tyc)+sum(((log(gridc(1:pc))-tavglc(tc,tyc))**2.0)*(cdist(tc,tyc,1:pc,sc,ic)/pitypes(tyc)))
		end do
	end do
	tstdlc(tc,tyc)=sqrt(tvarlc(tc,tyc))
end do

bvarlc(tc)=sum(  ( ( tavglc(tc,:)-avglc(tc) )**2.0 )*pitypes(:))
bstdlc(tc)=sqrt(bvarlc(tc))

open(unit=48,file='Decomp.txt')
do ttc=1,nt
	write(48,'(8f9.4)')  varlc(ttc), stdlc(ttc), bvarlc(ttc), bstdlc(ttc), tvarlc(ttc,:),tstdlc(ttc,:)
end do
close(48)


end subroutine condist

