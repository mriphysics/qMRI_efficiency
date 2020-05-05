function [eff,signals] = analysis_acq_set_PLANET(u, p, ...
    dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP,...
    signal_bSSFP)
%[eff,signals] = analysis_acq_set_PLANET(u, p, dmdT1_bSSFP, dmdT2_bSSFP, dmdM0_bSSFP, dmdP0_bSSFP, dmdB0_bSSFP, signal_bSSFP)
%   Calculates the efficiency and signals of a PLANET acquisition given the
%   acquisition and tissue parameters (u and p, respectively)

FIM = zeros(5);
aux_FIM = zeros(5);

nbSSFP = numel(u)-2;

signals = zeros(nbSSFP,1);

% CRLB and signal calculation
for ii=1:nbSSFP
    
    aux_u(1) = u(1);
    aux_u(2) = u(2);
    aux_u(3) = u(2+ii);
    
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

CRLB = inv(FIM); 

Tacq   = u(2)*nbSSFP*1e-3;               %acquisition time [s]
eff(1) = p(1) / sqrt(CRLB(1,1) * Tacq);  %T1 efficiency [s^-0.5]
eff(2) = p(2) / sqrt(CRLB(2,2) * Tacq);  %T2 efficiency [s^-0.5]

end

