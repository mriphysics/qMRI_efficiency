function [dmdO] = MEP_signal_derivatives(u, p, idx)
%[dmdO] = MEP_signal_derivatives(u, p, idx_p)
%   Calculates the derivatives of EPG states with respect to variables p
%   indicated in the indexes idx

stepsize = 1e-6;

% rows contain the states [F1plus; F0minus; F1minus], columns contain 
% different derivatives: dmdO(i,j) = d(Fi)/dOj
dmdO = zeros(3, numel(idx));

for ii=1:numel(idx)
    p_plus = p;
    p_plus(idx(ii)) = p_plus(idx(ii)) + stepsize;
    states_plus = MEP_signals(u,p_plus);
    
    p_minus = p;
    p_minus(idx(ii)) = p_minus(idx(ii)) - stepsize;
    states_minus = MEP_signals(u,p_minus);
    
    dmdO(:,ii) = (states_plus-states_minus)/(2*stepsize);
end

end

