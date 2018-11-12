subroutine borrowtr(x ,fv, n)
! this system deals with the case in which highest DC is binding

use params


implicit none
integer , intent(in) :: n
real(prec) , intent(in) ::x(n)
real(prec) , intent(out) :: fv(n)
real(prec)::consum,val,ass,vals(2)
integer:: inds(2)

! Consumption in complete markets in first state

ass=x(1)
call basefun (grida(tc+1,tyc,sc,ic,:),pa,ass,vals,inds)
val=vals(1)*vfun(tc+1,tyc,sc,ic,inds(1))+vals(2)*vfun(tc+1,tyc,sc,ic,inds(2))

fv(1)=val-vmaut(tc+1,tyc,sc,ic)

end subroutine borrowtr