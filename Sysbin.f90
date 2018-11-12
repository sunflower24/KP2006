subroutine sysbin(x ,fv, n)
! this system deals with the case in which highest DC is binding

use params


implicit none
integer , intent(in) :: n
real(prec) , intent(in) ::x(n)
real(prec) , intent(out) :: fv(n)
real(prec) :: ch,ah(ns,ni),npla(ns,ni),vals(2),erhs(n),elhs(n),valport(ns,ni)
integer inds(2),i,hc,hhc,helpc,helpcc


ah(1:ns,1:ni)=amin(tc,tyc,1:ns,1:ni)

do hc=1,hsc-1
	ah(hc,:)=x((hc-1)*ni+1:hc*ni)
end do
ah(hsc,1:hic)=x((hsc-1)*ni+1:n)

valport=0.0
do helpc=1,ns
	do helpcc=1,ni
		valport(helpc,helpcc)=(1.0/R(tc))*probs(tc,sc,helpc)*pi(tc,helpcc)*ah(helpc,helpcc) 
	end do
end do
 
ch=grida(tc,tyc,sc,ic,ac)+w(tc)*yat(tc,tyc,sc,ic) - sum(valport)


! Construct the Right Hand Side of the Euler Equations

do hc=1,hsc-1
	do hhc=1,ni
		call basefun (grida(tc,tyc,hc,hhc,:),pa,ah(hc,hhc),vals,inds)
		npla(hc,hhc)=vals(1)*lafun(tc,tyc,hc,hhc,inds(1))+vals(2)*lafun(tc,tyc,hc,hhc,inds(2))
		erhs((hc-1)*ni+hhc)=npla(hc,hhc)
	end do
end do
do hhc=1,hic
	call basefun (grida(tc,tyc,hsc,hhc,:),pa,ah(hsc,hhc),vals,inds)
	npla(hsc,hhc)=vals(1)*lafun(tc,tyc,hsc,hhc,inds(1))+vals(2)*lafun(tc,tyc,hsc,hhc,inds(2))
	erhs((hsc-1)*ni+hhc)=npla(hsc,hhc)
end do

! Construct Euler Equation

do i=1,n

	
	elhs(i)=(1.0-bta)*MU(ch)/(bta*R(tc))
	fv(i)=elhs(i)-erhs(i)

end do

end subroutine sysbin
