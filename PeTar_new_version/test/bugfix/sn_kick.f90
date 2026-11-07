program test_sn_kick
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer :: ns, i, j, kw, ecs, idum, bhflag, psflag, kmech, ecflag
  double precision :: mc, env, frac, total, expected, sigma, mxns
  double precision :: m1, m1n, m2, ecc, sep, jorb, v(4), mco
  double precision :: cores(7) = [2d0,3d0,5d0,6.5d0,8d0,10d0,12d0]
  double precision :: envelopes(6) = [0d0,1d-300,0.01d0,0.1d0,1d0,10d0]
  common /value3/ idum
  common /value4/ sigma,mxns,bhflag
  common /flags2/ psflag,kmech,ecflag
  do ns=3,4
    do i=1,size(cores)
      mc=cores(i)
      do j=1,size(envelopes)
        env=envelopes(j)
        call bounded_fallback(mc,env,ns,frac,total)
        if (.not.ieee_is_finite(frac+total)) stop 1
        if (frac<0d0.or.frac>1d0.or.total<0d0.or.total>env) stop 2
        if (abs(total-frac*env)>1d-14*max(env,1d-300)) stop 3
        if (mc==2d0.and.abs(total-min(0.2d0,env))>1d-15) stop 4
      end do
    end do
  end do
  ! Ordinary rapid/delayed cases retain the established prescriptions.
  call bounded_fallback(8d0,10d0,3,frac,total)
  expected=((-0.25d0+0.1275d0)*3d0+1d0)*10d0
  if(abs(total-expected)>1d-13) stop 5
  call bounded_fallback(5d0,10d0,4,frac,total)
  expected=((0.133d0-0.0093d0)*(5d0-11d0)+1d0)*10d0
  if(abs(total-expected)>1d-13) stop 6
  sigma=0d0; mxns=2.5d0; bhflag=2; psflag=0; kmech=1; ecflag=0
  do i=1,3
    idum=-12345; kw=13; ecs=0; m1=5d0; m1n=5d0; m2=3d0
    ecc=0.2d0; sep=1000d0; jorb=0d0; frac=0d0; total=0d0; mco=3d0
    if(i==2) sigma=265d0
    if(i==3) frac=1d0
    call kick(kw,m1,m1n,m2,ecc,sep,jorb,v,frac,total,mco,ecs)
    if(.not.all(ieee_is_finite(v))) stop 7
    if(.not.ieee_is_finite(ecc+sep+jorb)) stop 8
    if(i/=2.and.any(v/=0d0)) stop 9
    if(i==2.and.v(4)<=0d0) stop 10
    print *, 'kick case',i,v
  end do
  print *, 'fallback and zero-kick regressions passed'
end program
