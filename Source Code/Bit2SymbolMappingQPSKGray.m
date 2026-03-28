function d = Bit2SymbolMappingQPSKGray(A, m)
  %-------------------------------------------------------------------------
  % INPUT:  m    Massage 
  %         A    Amplitude of symbols
  % OUTPUT: d    Corresponding symbol vector
  % ----------------
  % (-1-j)/sqrt(2)    00
  % (-1+j)/sqrt(2)    01
  %  (1+j)/sqrt(2)    11
  %  (1-j)/sqrt(2)    10
  %-------------------------------------------------------------------------
d = (A/sqrt(2))*((2*m(1:2:end-1)-1) + (2*m(2:2:end)-1)*1j);
end 