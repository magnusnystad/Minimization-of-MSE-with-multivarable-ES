function y = smoothstep(x,k)
% Smootherstep function that can take in one or several arguments
% k scales x-axis of the function so that the step goes from x = 0 to x=k
x_mod = x/k;
y = zeros(size(x));
for i = 1:length(x)
    if x(i) <= 0
        y(i) = 0;
    elseif x(i) > 0 && x(i) < k
        y(i) = 3*x_mod(i)^2 - 2*x_mod(i)^3;
    else
        y(i) = 1;
    end
end
end

