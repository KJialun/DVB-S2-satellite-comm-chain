function [out] = RC(f, Ts, beta)
alpha = (1 - beta)/(2*Ts);
gamma = (1 + beta)/(2*Ts);
out = zeros(length(f),1);
for i=1:length(f)
    if abs(f(i)) <= alpha
        out(i) = Ts;
    elseif abs(f(i)) <= gamma && abs(f(i)) > alpha
        out(i) = (Ts/2)*(1 + cos((pi*Ts/beta)*(abs(f(i)) - alpha)));
    else
        out(i) = 0;
    end
end
end