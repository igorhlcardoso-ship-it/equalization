function z = AWGN(Delta, v, d)
  %-------------------------------------------------------------------------
  % INPUT:  Delta   Norm of impluse response h_e(t)
  %         v       Variance
  %         d       Symbol vector
  % OUTPUT: z       Vector of samples z(t0 + nTs)
  %-------------------------------------------------------------------------

  L = length(d);
  nl_real = sqrt(v)*randn(1,L);             
  nl_imag = sqrt(v)*randn(1,L);             
  nl = nl_real + 1i*nl_imag;

  z = ((Delta*Delta))*d + (Delta/sqrt(2))*nl;
end
