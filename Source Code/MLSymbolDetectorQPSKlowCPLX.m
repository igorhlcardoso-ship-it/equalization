function d_est = MLSymbolDetectorQPSKlowCPLX(A, Delta, z)
  %-------------------------------------------------------------------------
  % INPUT:  A       Amplitude of the constellation symbols
  %         Delta   Norm of impluse response h_e(t)
  %         z       Vector of samples z(t0 + nTs)
  % OUTPUT: d_est   The ML symbol estimated vector
  %-------------------------------------------------------------------------  
  
  d_est = (A/sqrt(2))*(sign(real(z)) + 1j*sign(imag(z)));

end

