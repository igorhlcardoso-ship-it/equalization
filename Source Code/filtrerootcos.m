function he = filtrerootcos(Beta, F, Omega)
  %-------------------------------------------------------------------------
  % USAGE:  Impulse response of root-raised-cosine filter from 0 to
  %         2*Omega*F (samples)
  % INPUT:  Beta    Roll-off factor of root-raised-cosine filter
  %         F       Number of samples / symbols
  % OUTPUT: he      Impulse response of root-raised-cosine filter
  %-------------------------------------------------------------------------
  % Author: Khoa HUYNH
  % Created: 09-01-2024
  %-------------------------------------------------------------------------
    tempf = ((1:2*Omega*F) - Omega*F)/F;
    he = (4*Beta/pi/sqrt(F)) * (cos((1+Beta)*pi*tempf) + ...
                      ((1-Beta)*pi/4/Beta)*sinc((1-Beta)*tempf)) ...
                      ./ (1 - (4*Beta*tempf).^2);
    he(Omega*F) == 1;
end
