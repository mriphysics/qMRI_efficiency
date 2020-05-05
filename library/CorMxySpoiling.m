function [C50] = CorMxySpoiling(FA,TR)
%[C50] = CorMxySpoiling(FA,TR)
%   Correction for incomplete RF spoiling according to Baudrexel et al 2018
%   (DOI:10.1002/mrm21120)

P50 = [9.639e-1     4.989e-3        -1.254e-4       -3.18e-6        1.527e-7        -1.462e-9;
       5.88e-3     -1.056e-3         4.801e-5       -8.549e-7       5.382e-9         0;
       4.143e-4    -4.92e-6         -1.56e-7         2.282e-9       0                0;
      -1.5059e-5    2.334e-7        -1.189e-9        0              0                0;
       9.449e-8    -1.025e-9         0               0              0                0;
      -4.255e-10    0                0               0              0                0];
   
TR_poly = TR*ones(6,1);
FA_poly = FA*ones(1,6);
for ii=1:6
    TR_poly(ii) = (TR_poly(ii))^(ii-1);
    FA_poly(ii) = (FA_poly(ii))^(ii-1);   
end

C50 = FA_poly * P50 * TR_poly;

