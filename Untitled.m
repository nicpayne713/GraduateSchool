f5 = @(x) exp(2.*x)./((1+x.^2).^2);
q = quad(f5, 0,3)
ql = quadl(f5,0,3)
qgk = quadgk(f5,0,3)
int = integral(f5, 0,3)