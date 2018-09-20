function c = t_project_coefficients_ab ( n, f, a, b )
%  Parameters:
%
%    Input, integer N, the highest order polynomial to compute.
%
%    Input, function handle F, of the form
%      function v = f ( x )
%
%    Input, real A, B, the interval of definition.
%
%    Output, real C(N+1,1), the projection coefficients of f(x) onto
%    T(0,x) through T(n,x), based on [A,B].
%
  for k = 1 : n + 1

    t = cos ( pi * ( k - 0.5 ) / ( n + 1 ) ); %recursive definition of Chebyshev Polynomial

    y = ( ( 1.0 + t ) * b   ...
        + ( 1.0 - t ) * a ) ...
        /   2.0; %maps each point in [a,b] to [-1,1]

    d(k) = f ( y );

  end

  fac = 2.0 / ( n + 1 );
  for j = 1 : ( n + 1 )
    sum = 0.0;
    for k = 1 : ( n + 1 )
      sum = sum + d(k) * cos ( ( pi * ( j - 1 ) ) * ( ( k - 0.5 ) / ( n + 1 ) ) );
    end
    c(j) = fac * sum;
  end

  c(1) = c(1) / 2.0;

  return
end