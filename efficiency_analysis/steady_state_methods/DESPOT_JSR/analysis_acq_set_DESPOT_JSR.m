function [eff, signals] = analysis_acq_set_DESPOT_JSR(u, p, nbSSFP, ...
    dmdT1_SPGR,  dmdT2_SPGR,  dmdM0_SPGR,  dmdP0_SPGR,  dmdB0_SPGR, ...
    dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP,...
    signal_SPGR, signal_bSSFP)
%[eff, signals] = analysis_acq_set_DESPOT_JSR(u, p, nbSSFP, dmdT1_SPGR, dmdT2_SPGR, dmdM0_SPGR, dmdP0_SPGR, dmdB0_SPGR, dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP, signal_SPGR, signal_bSSFP)
%   Calculates the efficiency and signals of a DESPOT/JSR acquisition given 
%   the acquisition and tissue parameters (u and p, respectively)

FIM = zeros(5);
aux_FIM = zeros(5);

nSS = numel(u)/3; %number of steady-states
nSPGR = nSS - nbSSFP;

signals = zeros(nSS,1);

% CRLB and signal calculation
for ii=1:nSS
    if ii<=nSPGR

        aux_u(1) = u(ii);
        aux_u(2) = u(ii + nSS);
        aux_u(3) = u(ii + 2*nSS);
        
        signals(ii) = signal_SPGR(aux_u,p);
     
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
        
        signals(ii) = signal_bSSFP(aux_u,p);
        
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

CRLB = inv(FIM);

Tacq   = sum(u(nSS+1:2*nSS))*1e-3;       %acquisition time [s]
eff(1) = p(1) / sqrt(CRLB(1,1) * Tacq);  %T1 efficiency [s^-1/2]
eff(2) = p(2) / sqrt(CRLB(2,2) * Tacq);  %T2 efficiency [s^-1/2]

end

