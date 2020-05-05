function [A] = gen_comb_steady_states(max_nSS,min_nSS)
%[A] = gen_comb_steady_states(max_nSS,min_nSS)
%   Generates combinations of steady-states with a minimum of min_nSS
%   steady-states and a maximum of max_nSS steady-states

% generates all combinations of 2 elements
aux_comb = combnk(0:max_nSS,2); 
% remove combinations with more than max_nSS steady-states
aux_comb(sum(aux_comb,2)>max_nSS,:) = []; 
aux_comb(sum(aux_comb,2)<min_nSS,:) = []; 

% find all permutations of the selected combinations
aux_2_comb = [];
for ii=1:size(aux_comb,1)
    a = perms(aux_comb(ii,:));
    aux_2_comb = [aux_2_comb; a];
end

% remove possible repeated combinations of steady-states
A = unique(aux_2_comb,'rows');

% add combinations that have the same number of SPGR and bSSFP
for ii=ceil(min_nSS/2):floor(max_nSS/2)
    A = [A; [ii ii]];
end


end

