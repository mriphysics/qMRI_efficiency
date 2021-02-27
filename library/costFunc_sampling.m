function [hnd_sampling] = costFunc_sampling(hnd_CostFunc,P,varargin)
%Samples the cost function with equal weights (further extension would be 
%receiving weighting matrix W) for different parameters (e.g. several T1 
%and T2 combinations)
%The final cost function is then the sum of all samples
%   handle_costFunc is a handle to the cost function that receives two 
%inputs:
%1 - the design variable sused for the optimization
%2 - the parameter values at which the cost function should be sampled
%   P is a NxL matrix where N is the number of sample sof the cost function
%and L is the number of different parameters that the second input of
%hand_costFunc should receive every time
%   handle_sampling is a handle to a function which only input is the same
%   as input 1 for handle_costFunc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   David Leitao; 10/09/2018 
%   david.leitao@kcl.ac.uk


%Check optional inputs
%...

numb_samples = size(P,1);

hnd_sampling =@(x) hnd_CostFunc(x,P(1,:));

ii = 1;
while ii<numb_samples 
    
    ii = ii + 1;
    
    hnd_sampling =@(x) hnd_CostFunc(x,P(ii,:)) + hnd_sampling(x);
    
%     Optional max or min
%     handle_sampling =@(x) max([handle_costFunc(x,P(ii,:)), handle_sampling(x)]);
    
end

end

