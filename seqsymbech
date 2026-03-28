function dech = seqsymbech(d, F)
  %-------------------------------------------------------------------------
  % USAGE:  Sample continuous signal d <=> insert F-1 zero element into
  %         vector d
  % INPUT:  d       Symbol sequences after mapping
  %         F       Number of samples / symbols
  % OUTPUT: dech    Corresponding sample vector
  %-------------------------------------------------------------------------
   dech = zeros(1, F * length(d));
   dech(1:F:end-F+1) = d;
end
