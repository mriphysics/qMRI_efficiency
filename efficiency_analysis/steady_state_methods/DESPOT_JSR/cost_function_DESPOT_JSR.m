function [f] = cost_function_DESPOT_JSR(u, p, nbSSFP, ...
    dmdT1_SPGR,  dmdT2_SPGR,  dmdM0_SPGR,  dmdP0_SPGR,  dmdB0_SPGR, ...
    dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP)
%[f] = cost_function_DESPOT_JSR(u, p, nbSSFP, dmdT1_SPGR, dmdT2_SPGR, dmdM0_SPGR, dmdP0_SPGR, dmdB0_SPGR, dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP)
%   Calculates cost function of DESPOT/JSR for a single parameter vector p

FIM = zeros(5);
aux_FIM = zeros(5);

nSS = numel(u)/3; %number of steady-states
nSPGR = nSS - nbSSFP;

% CRLB calculation
for ii=1:nSS
    if ii<=nSPGR
        
        aux_u(1) = u(ii);
        aux_u(2) = u(ii + nSS);
        aux_u(3) = u(ii + 2*nSS);
     
        aux_FIM(1,1) = real(dmdT1_SPGR(aux_u,p)) * real(dmdT1_SPGR(aux_u,p));
        aux_FIM(1,2) = real(dmdT1_SPGR(aux_u,p)) * real(dmdT2_SPGR(aux_u,p)); 
        aux_FIM(1,3) = real(dmdT1_SPGR(aux_u,p)) * real(dmdM0_SPGR(aux_u,p));
        aux_FIM(1,4) = real(dmdT1_SPGR(aux_u,p)) * real(dmdP0_SPGR(aux_u,p));
        aux_FIM(1,5) = real(dmdT1_SPGR(aux_u,p)) * real(dmdB0_SPGR(aux_u,p));
        aux_FIM(2,2) = real(dmdT2_SPGR(aux_u,p)) * real(dmdT2_SPGR(aux_u,p));
        aux_FIM(2,3) = real(dmdT2_SPGR(aux_u,p)) * real(dmdM0_SPGR(aux_u,p));
        aux_FIM(2,4) = real(dmdT2_SPGR(aux_u,p)) * real(dmdP0_SPGR(aux_u,p));
        aux_FIM(2,5) = real(dmdT2_SPGR(aux_u,p)) * real(dmdB0_SPGR(aux_u,p));
        aux_FIM(3,3) = real(dmdM0_SPGR(aux_u,p)) * real(dmdM0_SPGR(aux_u,p));
        aux_FIM(3,4) = real(dmdM0_SPGR(aux_u,p)) * real(dmdP0_SPGR(aux_u,p));
        aux_FIM(3,5) = real(dmdM0_SPGR(aux_u,p)) * real(dmdB0_SPGR(aux_u,p));
        aux_FIM(4,4) = real(dmdP0_SPGR(aux_u,p)) * real(dmdP0_SPGR(aux_u,p));
        aux_FIM(4,5) = real(dmdP0_SPGR(aux_u,p)) * real(dmdB0_SPGR(aux_u,p));
        aux_FIM(5,5) = real(dmdB0_SPGR(aux_u,p)) * real(dmdB0_SPGR(aux_u,p));
        %symmetric terms
        aux_FIM(2,1) = aux_FIM(1,2);
        aux_FIM(3,1) = aux_FIM(1,3);
        aux_FIM(4,1) = aux_FIM(1,4);
        aux_FIM(5,1) = aux_FIM(1,5);
        aux_FIM(3,2) = aux_FIM(2,3);
        aux_FIM(4,2) = aux_FIM(2,4);
        aux_FIM(5,2) = aux_FIM(2,5);
        aux_FIM(4,3) = aux_FIM(3,4);
        aux_FIM(5,3) = aux_FIM(3,5);
        aux_FIM(5,4) = aux_FIM(4,5);
        
        FIM = FIM + aux_FIM;
        
        aux_FIM(1,1) = imag(dmdT1_SPGR(aux_u,p)) * imag(dmdT1_SPGR(aux_u,p));
        aux_FIM(1,2) = imag(dmdT1_SPGR(aux_u,p)) * imag(dmdT2_SPGR(aux_u,p)); 
        aux_FIM(1,3) = imag(dmdT1_SPGR(aux_u,p)) * imag(dmdM0_SPGR(aux_u,p));
        aux_FIM(1,4) = imag(dmdT1_SPGR(aux_u,p)) * imag(dmdP0_SPGR(aux_u,p));
        aux_FIM(1,5) = imag(dmdT1_SPGR(aux_u,p)) * imag(dmdB0_SPGR(aux_u,p));
        aux_FIM(2,2) = imag(dmdT2_SPGR(aux_u,p)) * imag(dmdT2_SPGR(aux_u,p));
        aux_FIM(2,3) = imag(dmdT2_SPGR(aux_u,p)) * imag(dmdM0_SPGR(aux_u,p));
        aux_FIM(2,4) = imag(dmdT2_SPGR(aux_u,p)) * imag(dmdP0_SPGR(aux_u,p));
        aux_FIM(2,5) = imag(dmdT2_SPGR(aux_u,p)) * imag(dmdB0_SPGR(aux_u,p));
        aux_FIM(3,3) = imag(dmdM0_SPGR(aux_u,p)) * imag(dmdM0_SPGR(aux_u,p));
        aux_FIM(3,4) = imag(dmdM0_SPGR(aux_u,p)) * imag(dmdP0_SPGR(aux_u,p));
        aux_FIM(3,5) = imag(dmdM0_SPGR(aux_u,p)) * imag(dmdB0_SPGR(aux_u,p));
        aux_FIM(4,4) = imag(dmdP0_SPGR(aux_u,p)) * imag(dmdP0_SPGR(aux_u,p));
        aux_FIM(4,5) = imag(dmdP0_SPGR(aux_u,p)) * imag(dmdB0_SPGR(aux_u,p));
        aux_FIM(5,5) = imag(dmdB0_SPGR(aux_u,p)) * imag(dmdB0_SPGR(aux_u,p));
        %symmetric terms
        aux_FIM(2,1) = aux_FIM(1,2);
        aux_FIM(3,1) = aux_FIM(1,3);
        aux_FIM(4,1) = aux_FIM(1,4);
        aux_FIM(5,1) = aux_FIM(1,5);
        aux_FIM(3,2) = aux_FIM(2,3);
        aux_FIM(4,2) = aux_FIM(2,4);
        aux_FIM(5,2) = aux_FIM(2,5);
        aux_FIM(4,3) = aux_FIM(3,4);
        aux_FIM(5,3) = aux_FIM(3,5);
        aux_FIM(5,4) = aux_FIM(4,5);
        
        FIM = FIM + aux_FIM;
        
    else %bSSFP sequence
        
        aux_u(1) = u(ii);
        aux_u(2) = u(ii + nSS);
        aux_u(3) = u(ii + 2*nSS);
        
        aux_FIM(1,1) = real(dmdT1_bSSFP(aux_u,p)) * real(dmdT1_bSSFP(aux_u,p));
        aux_FIM(1,2) = real(dmdT1_bSSFP(aux_u,p)) * real(dmdT2_bSSFP(aux_u,p)); 
        aux_FIM(1,3) = real(dmdT1_bSSFP(aux_u,p)) * real(dmdM0_bSSFP(aux_u,p));
        aux_FIM(1,4) = real(dmdT1_bSSFP(aux_u,p)) * real(dmdP0_bSSFP(aux_u,p));
        aux_FIM(1,5) = real(dmdT1_bSSFP(aux_u,p)) * real(dmdB0_bSSFP(aux_u,p));
        aux_FIM(2,2) = real(dmdT2_bSSFP(aux_u,p)) * real(dmdT2_bSSFP(aux_u,p));
        aux_FIM(2,3) = real(dmdT2_bSSFP(aux_u,p)) * real(dmdM0_bSSFP(aux_u,p));
        aux_FIM(2,4) = real(dmdT2_bSSFP(aux_u,p)) * real(dmdP0_bSSFP(aux_u,p));
        aux_FIM(2,5) = real(dmdT2_bSSFP(aux_u,p)) * real(dmdB0_bSSFP(aux_u,p));
        aux_FIM(3,3) = real(dmdM0_bSSFP(aux_u,p)) * real(dmdM0_bSSFP(aux_u,p));
        aux_FIM(3,4) = real(dmdM0_bSSFP(aux_u,p)) * real(dmdP0_bSSFP(aux_u,p));
        aux_FIM(3,5) = real(dmdM0_bSSFP(aux_u,p)) * real(dmdB0_bSSFP(aux_u,p));
        aux_FIM(4,4) = real(dmdP0_bSSFP(aux_u,p)) * real(dmdP0_bSSFP(aux_u,p));
        aux_FIM(4,5) = real(dmdP0_bSSFP(aux_u,p)) * real(dmdB0_bSSFP(aux_u,p));
        aux_FIM(5,5) = real(dmdB0_bSSFP(aux_u,p)) * real(dmdB0_bSSFP(aux_u,p));
        %symmetric terms
        aux_FIM(2,1) = aux_FIM(1,2);
        aux_FIM(3,1) = aux_FIM(1,3);
        aux_FIM(4,1) = aux_FIM(1,4);
        aux_FIM(5,1) = aux_FIM(1,5);
        aux_FIM(3,2) = aux_FIM(2,3);
        aux_FIM(4,2) = aux_FIM(2,4);
        aux_FIM(5,2) = aux_FIM(2,5);
        aux_FIM(4,3) = aux_FIM(3,4);
        aux_FIM(5,3) = aux_FIM(3,5);
        aux_FIM(5,4) = aux_FIM(4,5);
        
        FIM = FIM + aux_FIM;
        
        aux_FIM(1,1) = imag(dmdT1_bSSFP(aux_u,p)) * imag(dmdT1_bSSFP(aux_u,p));
        aux_FIM(1,2) = imag(dmdT1_bSSFP(aux_u,p)) * imag(dmdT2_bSSFP(aux_u,p)); 
        aux_FIM(1,3) = imag(dmdT1_bSSFP(aux_u,p)) * imag(dmdM0_bSSFP(aux_u,p));
        aux_FIM(1,4) = imag(dmdT1_bSSFP(aux_u,p)) * imag(dmdP0_bSSFP(aux_u,p));
        aux_FIM(1,5) = imag(dmdT1_bSSFP(aux_u,p)) * imag(dmdB0_bSSFP(aux_u,p));
        aux_FIM(2,2) = imag(dmdT2_bSSFP(aux_u,p)) * imag(dmdT2_bSSFP(aux_u,p));
        aux_FIM(2,3) = imag(dmdT2_bSSFP(aux_u,p)) * imag(dmdM0_bSSFP(aux_u,p));
        aux_FIM(2,4) = imag(dmdT2_bSSFP(aux_u,p)) * imag(dmdP0_bSSFP(aux_u,p));
        aux_FIM(2,5) = imag(dmdT2_bSSFP(aux_u,p)) * imag(dmdB0_bSSFP(aux_u,p));
        aux_FIM(3,3) = imag(dmdM0_bSSFP(aux_u,p)) * imag(dmdM0_bSSFP(aux_u,p));
        aux_FIM(3,4) = imag(dmdM0_bSSFP(aux_u,p)) * imag(dmdP0_bSSFP(aux_u,p));
        aux_FIM(3,5) = imag(dmdM0_bSSFP(aux_u,p)) * imag(dmdB0_bSSFP(aux_u,p));
        aux_FIM(4,4) = imag(dmdP0_bSSFP(aux_u,p)) * imag(dmdP0_bSSFP(aux_u,p));
        aux_FIM(4,5) = imag(dmdP0_bSSFP(aux_u,p)) * imag(dmdB0_bSSFP(aux_u,p));
        aux_FIM(5,5) = imag(dmdB0_bSSFP(aux_u,p)) * imag(dmdB0_bSSFP(aux_u,p));
        %symmetric terms
        aux_FIM(2,1) = aux_FIM(1,2);
        aux_FIM(3,1) = aux_FIM(1,3);
        aux_FIM(4,1) = aux_FIM(1,4);
        aux_FIM(5,1) = aux_FIM(1,5);
        aux_FIM(3,2) = aux_FIM(2,3);
        aux_FIM(4,2) = aux_FIM(2,4);
        aux_FIM(5,2) = aux_FIM(2,5);
        aux_FIM(4,3) = aux_FIM(3,4);
        aux_FIM(5,3) = aux_FIM(3,5);
        aux_FIM(5,4) = aux_FIM(4,5);
        
        FIM = FIM + aux_FIM;
        
    end
end

if cond(FIM)>1e15
    % if condition number is too high inv() might return unrealistic 
    % negative/complex numbers
    f = Inf;
else
    CRLB = inv(FIM);

    Tacq = sum(u(nSS+1:2*nSS))*1e-3;          %acquisition time [s]
    sq_eff_T1 = p(1)^2 / (CRLB(1,1) * Tacq);  %squared T1 efficiency [s^-1]
    sq_eff_T2 = p(2)^2 / (CRLB(2,2) * Tacq);  %squared T2 efficiency [s^-1]

    % cost function is the sum of squares of 1/efficiency
    f = 1/sq_eff_T1 + 1/sq_eff_T2;
end

end

