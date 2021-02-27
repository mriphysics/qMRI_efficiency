function [signals] = MEP_signals(u,p)
%[signals] = MEP_signals(u,p)
%   Calculates Multiple Echo Pathways steady-state signals of spoiled GRE 
%   with no RF spoiling for F1plus, F0minus and F1minus states

%u = [FA TR, TE0, dTE1, dTE2]
%p = [T1 T2 M0]

N = 1000; 
FA = u(1)*ones(N,1);
RF = zeros(N,1);

TE0 = u(3);
TE1 = u(3)+u(4);
TE2 = u(3)+u(4)+u(5);

E2_TEp1 = exp(-TE0/p(2)); %TE for F1plus state
E2_TEm0 = exp(-TE1/p(2)); %TE for F0minus state
E2_TEm1 = exp(-TE2/p(2)); %TE for F1minus state

%                              [N, FA, RF, TR,   T1,   T2,   M0]
[~, states] = cppEPG_GRE_MEP(N, FA, RF, u(2), p(1), p(2), p(3));

% rows contain the states [F1plus; F0minus; F1minus] and each column is a
% different (dephasing) k order: {0,1,...,24}
F1plus  = states(1,2) * E2_TEp1;
F0minus = states(2,1) * E2_TEm0;
F1minus = states(2,2) * E2_TEm1;

signals = [F1plus; 
           F0minus; 
           F1minus];

end

