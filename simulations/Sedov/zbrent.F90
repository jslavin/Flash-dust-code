real function zbrent(func,x1,x2,tol)
    !! root finding routine based on Brent's method (via Numerical Recipes)
    !! translated from Fortran 77 version
    implicit none
    external func
    real func
    real, intent(IN) :: x1,x2,tol
    real a,b,c,fa,fb,fc,eps,d,e,tol1,xm,p,q,r,s
    integer iter,itmax
    parameter (itmax=100,eps=3.E-8)
    a = x1
    b = x2
    fa = func(a)
    fb = func(b)
    if(fb*fa .gt. 0.) then
        write(*,*) 'Root must be bracketed for ZBRENT.'
        zbrent = 0.
        return
    endif
    fc = fb
    do iter=1,itmax
        if(fb*fc .gt. 0.) then
            c = a
            fc = fa
            d = b - a
            e = d
        endif
        if(ABS(fc) .lt. ABS(fb)) then
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        endif
        tol1 = 2.*eps*ABS(b) + 0.5*tol
        xm = .5*(c - b)
        if((ABS(xm) .le. tol1) .or. (fb .eq. 0.)) then
            zbrent = b
            fb = func(b)
            return
        endif
        if((ABS(e) .ge. tol1) .and. (ABS(fa) .gt. ABS(fb))) then
            s = fb/fa
            if(a .eq. c) then
                p = 2.*xm*s
                q = 1. - s
            else
                q = fa/fc
                r = fb/fc
                p = s*(2.*xm*q*(q - r) - (b - a)*(r - 1.))
                q = (q - 1.)*(r - 1.)*(s - 1.)
            endif
            if(p .gt. 0.) q = -q
            p = ABS(p)
            if(2.*p  .lt. MIN(3.*xm*q - ABS(tol1*q),ABS(e*q))) then
                e = d
                d = p/q
            else
                d = xm
                e = d
            endif
        else
            d = xm
            e = d
        endif
        a = b
        fa = fb
        if(ABS(d) .gt. tol1) then
            b = b + d
        else
            b = b + SIGN(tol1,xm)
        endif
        fb = func(b)
    end do
    write(*,*) 'ZBRENT exceeding maximum iterations.'
    zbrent = b
    return
end function zbrent
