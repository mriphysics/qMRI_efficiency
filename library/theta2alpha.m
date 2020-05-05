function [alpha] = theta2alpha(theta)
%[alpha] = theta2alpha(theta)

alpha = zeros(size(theta));

alpha(1) = theta(1);
alpha(2:end) = theta(2:end) + theta(1:end-1);
end

