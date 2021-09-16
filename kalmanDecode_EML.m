%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


function x_est = kalmanDecode_EML(z,x_prev,A,H,K)
  %#codegen 
  
  %kalmanDecode_eml:

  x_est = A*x_prev + K *(z -H*A*x_prev);
  
end