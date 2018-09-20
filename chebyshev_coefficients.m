function c = chebyshev_coefficients ( a, b, n, f )
%  Parameters:
%
%    Input, real A, B, the domain of definition.
%
%    Input, integer N, the order of the interpolant.
%
%    Input, real F ( X ), a function handle.
%
%    Output, real C(N), the Chebyshev coefficients.
%
   angle = [( 1 : 2 : ( 2 * n - 1 ) ) * pi / ( 2 * n )]';
  x = cos ( angle );
  x = 0.5 * ( a + b ) + x * 0.5 * ( b - a );
  fx = f ( x );

  c = zeros ( n, 1 );
  for j = 1 : n
    c(j) = 0.0;
    for k = 1 : n
      c(j) = c(j) + fx(k) * cos ( pi * ( j - 1 ) * ( 2 * k - 1 ) / 2 / n );
    end
  end

  c(1:n) = 2.0 * c(1:n) / n;

  return
end