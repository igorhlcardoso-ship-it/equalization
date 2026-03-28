function n = ErrosBit(m, m_est)
  %-----------------------------------------------------------------------------
  % INPUT:  m       Massage, vector of k = (kappa x L) bits to be sent
  %         m_est   The estimated vector of m
  % OUTPUT: n       The number of errous bits
  %-----------------------------------------------------------------------------
   k = length(m);           %Number of bits k
   n = 0;
   for i=1:k
     if abs(m_est(i) - m(i))==1
       n = n + 1;
     end
   end
end