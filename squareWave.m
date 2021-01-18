function y = squareWave(t,t0,period)
% Function that takes in current time, t, "start of wave time", t0, and
% period of the square wave (all in same unit). func gives current value of
% a unit amplitude square wave.

y = zeros(1,length(t));
    for i = 1:length(t)
        y(i) = 2*(2*floor((t(i)-t0)/period)-floor(2*(t(i)-t0)/period))+1;
    end
    
end

