function f=check_positive(u,func);

if any(u<-1e-6)
  save udump1 u
  disp(['NEGATIVE FOUND1! ' num2str(min(u))])
  f=func(u);
elseif any(u<-1e-4)
  save udump2 u
  disp(['NEGATIVE FOUND2! ' num2str(min(u))])
  f=func(u);
elseif any(u<-1e-2)
  save udump3 u
  disp(['NEGATIVE FOUND3! ' num2str(min(u))])
  f=func(u);
elseif any(u<0)
  save udump4 u
  disp(['NEGATIVE FOUND4! ' num2str(min(u))])
  f=func(u);
  save udump4 u
else
  save udump5 u
  f=func(u);
end

if any(~isfinite(f))
  save unan u f
  error('NaN detected!')
end
