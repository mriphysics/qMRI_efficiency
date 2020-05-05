function [f] = cost_function_DESPOT1(u, p, nSPGR,...
    dmdT1_SPGR,  dmdM0_SPGR)
%[f] = cost_function_DESPOT1(u, p, nSPGR, dmdT1_SPGR,  dmdM0_SPGR)
%   Calculates the cost function for DESPOT1 optimisation (inverse of T1
%   efficiency)

FIM = zeros(2);
aux_FIM = zeros(2);

for i1=1:nSPGR
    aux_u(1) = u(i1);
    aux_u(2) = u(i1 + nSPGR);

    aux_FIM(1,1) = real(dmdT1_SPGR(aux_u,p)) * real(dmdT1_SPGR(aux_u,p));
    aux_FIM(1,2) = real(dmdT1_SPGR(aux_u,p)) * real(dmdM0_SPGR(aux_u,p)); 
    aux_FIM(2,2) = real(dmdM0_SPGR(aux_u,p)) * real(dmdM0_SPGR(aux_u,p));
    %symmetric terms
    aux_FIM(2,1) = aux_FIM(1,2);
    
    FIM = FIM + aux_FIM;
    

    aux_FIM(1,1) = imag(dmdT1_SPGR(aux_u,p)) * imag(dmdT1_SPGR(aux_u,p));
    aux_FIM(1,2) = imag(dmdT1_SPGR(aux_u,p)) * imag(dmdM0_SPGR(aux_u,p)); 
    aux_FIM(2,2) = imag(dmdM0_SPGR(aux_u,p)) * imag(dmdM0_SPGR(aux_u,p));
    %symmetric terms
    aux_FIM(2,1) = aux_FIM(1,2);
    
    FIM = FIM + aux_FIM;
end


if cond(FIM)>1e15
    % if condition number is too high inv() might return unrealistic 
    % negative/complex numbers
    f = Inf;
else
    CRLB = inv(FIM);

    Tacq = sum(u(nSPGR+1:2*nSPGR))*1e-3;     %acquisition time [s]
    sq_eff_T1 = p(1)^2 / (CRLB(1,1) * Tacq); %squared T1 efficiency [s^-1]

    % cost function is 1/efficiency
    f = 1/sq_eff_T1;
end

end

