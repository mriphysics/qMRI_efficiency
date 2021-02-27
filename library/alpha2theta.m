function [theta] = alpha2theta(alpha)
%[theta] = alpha2theta(alpha)

theta = zeros(size(alpha));

theta(1) = alpha(1);
theta(2:end) = alpha(2:end)/2;

end

