function [eff,signals] = analysis_acq_set_DESS(u, p)
%[eff,signals] = analysis_acq_set_DESS(u, p)
%   Calculates the efficiency and signals of a DESS acquisition given the
%   acquisition and tissue parameters (u and p, respectively)

FIM = zeros(3);
aux_FIM = zeros(3);

nDESS   = length(u)/4;
nEchoes = 2;

signals = zeros(3,nDESS); %initially save all pathways, then discard those not measured

% CRLB and signal calculation
for ii=1:nDESS
    
    aux_u    = zeros(5,1);
    aux_u(1) = u(ii);
    aux_u(2) = u(ii+nDESS);
    aux_u(3) = 0;             %echo corresponding to F1plus is not used
    aux_u(4) = u(ii+2*nDESS); %TE0 corresponds to F0minus
    aux_u(5) = u(ii+3*nDESS); %dTE1 corresponds to F1minus
    
    signals(:,ii) = MEP_signals(aux_u,p);
    
    %dmdO(i,j) = d(Fi)/dOj, where Fi = [F1plus; F0minus; F1minus]
    dmdO = MEP_signal_derivatives(aux_u, p, 1:3);
    %DESS measures F0minus and F1minus:
    dmdT1 = dmdO(2:3,1); dmdT2 = dmdO(2:3,2); dmdM0 = dmdO(2:3,3); 
    
    for i2=1:nEchoes
        aux_FIM(1,1) = real(dmdT1(i2)) * real(dmdT1(i2));
        aux_FIM(1,2) = real(dmdT1(i2)) * real(dmdT2(i2)); 
        aux_FIM(1,3) = real(dmdT1(i2)) * real(dmdM0(i2));
        aux_FIM(2,2) = real(dmdT2(i2)) * real(dmdT2(i2));
        aux_FIM(2,3) = real(dmdT2(i2)) * real(dmdM0(i2));
        aux_FIM(3,3) = real(dmdM0(i2)) * real(dmdM0(i2));
        %symmetric terms
        aux_FIM(2,1) = aux_FIM(1,2);
        aux_FIM(3,1) = aux_FIM(1,3);
        aux_FIM(3,2) = aux_FIM(2,3);

        FIM = FIM + aux_FIM;
        
        aux_FIM(1,1) = imag(dmdT1(i2)) * imag(dmdT1(i2));
        aux_FIM(1,2) = imag(dmdT1(i2)) * imag(dmdT2(i2)); 
        aux_FIM(1,3) = imag(dmdT1(i2)) * imag(dmdM0(i2));
        aux_FIM(2,2) = imag(dmdT2(i2)) * imag(dmdT2(i2));
        aux_FIM(2,3) = imag(dmdT2(i2)) * imag(dmdM0(i2));
        aux_FIM(3,3) = imag(dmdM0(i2)) * imag(dmdM0(i2));
        %symmetric terms
        aux_FIM(2,1) = aux_FIM(1,2);
        aux_FIM(3,1) = aux_FIM(1,3);
        aux_FIM(3,2) = aux_FIM(2,3);

        FIM = FIM + aux_FIM;
    end
end

signals = signals(2:3,:); %discard pathway not measured (F1plus)

CRLB = inv(FIM);

Tacq = sum(u(nDESS+1:2*nDESS))*1e-3;     %acquisition time [s]
eff(1) = p(1) / sqrt(CRLB(1,1) * Tacq);  %T1 efficiency [s^-0.5]
eff(2) = p(2) / sqrt(CRLB(2,2) * Tacq);  %T2 efficiency [s^-0.5]

end

