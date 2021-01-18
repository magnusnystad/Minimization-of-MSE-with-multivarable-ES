function y = sat(s,b)
% Function that takes in an argument (s) and a ± limit (b) that determines
% the "type" of output (= sign(s) for s > b, s/b for s < b). Essentially a
% modified sign function that switches between 1 and -1 more gradually
    if abs(s) >= b
        y = sign(s);
    else
        y = s/b;
    end
    
end

