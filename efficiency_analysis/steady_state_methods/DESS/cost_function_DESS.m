function [f] = cost_function_DESS(u, p)
%[f] = cost_function_DESS(u, p)
%   Calculates cost function of DESS for a single parameter vector p

FIM = zeros(3);
aux_FIM = zeros(3);

nDESS   = numel(u)/4;
nEchoes = 2;

% CRLB calculation
for i1=1:nDESS
    
    aux_u    = zeros(5,1);
    aux_u(1) = u(i1);
    aux_u(2) = u(i1+nDESS);
    aux_u(3) = 0;             %echo corresponding to F1plus is not used 
    aux_u(4) = u(i1+2*nDESS); %TE0 corresponds to F0minus
    aux_u(5) = u(i1+3*nDESS); %dTE1 corresponds to F1minus
    
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

if cond(FIM)>1e15
    % if condition number is too high inv() might return unrealistic 
    % negative/complex numbers
    f = Inf;
else
    CRLB = inv(FIM);

    Tacq = sum(u(nDESS+1:2*nDESS))*1e-3;      %acquisition time [s]
    sq_eff_T1 = p(1)^2 / (CRLB(1,1) * Tacq);  %squared T1 efficiency [s^-1]
    sq_eff_T2 = p(2)^2 / (CRLB(2,2) * Tacq);  %squared T2 efficiency [s^-1]

    % cost function is the sum of squares of 1/efficiency
    f = 1/sq_eff_T1 + 1/sq_eff_T2;
end

end

