subroutine borrow(x ,fv, n)

use params


implicit none
integer , intent(in) :: n
real(prec) , intent(in) ::x(n)
real(prec) , intent(out) :: fv(n)
real(prec)::consum,val,ass,vals(2)
integer:: inds(2)

ass=x(1)
call basefun (grida(tc,tyc,sc,ic,:),pa,ass,vals,inds)
val=vals(1)*vfun(tc,tyc,sc,ic,inds(1))+vals(2)*vfun(tc,tyc,sc,ic,inds(2))

fv(1)=val-vmaut(tc,tyc,sc,ic)

end subroutine borrow
