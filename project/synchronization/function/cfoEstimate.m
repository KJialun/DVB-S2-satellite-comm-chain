function [n, df] = cfoEstimate(rx, p, T, K)
%{
    input
    rx: received symbols
    p:  pilot symbols
    T:  symbol duration
    K:  K-factor

    output
    n:  arrival time estimate
    df: cfo estimate
%}

N = length(p);  %% pilot symbols
L = length(rx); 
D = zeros(K, L-N+1); %% L-N +1 = length of data , 

for k = 1 : K
    for n = 1 : L - N 
        tmp = 0;
        for l = k : N-1
            tmp1 = conj(rx(n+l+1)) * p(l+1);
            
            tmp2 = conj(rx(n+l-k+1)) * p(l-k+1);
            
            tmp = tmp + tmp1 * conj(tmp2);
        end
        D(k,n) = tmp; %% different k size
        
    end 
    D(k,:) = D(k,:)/(N-k);
end

% Time of arrival estimate
tmp = sum(abs(D),1);
[~, n] = max(tmp);

% CFO estimate

k = (1:K).' ;
df = sum(angle(D(k,n))./(2*pi*k*T), 1);
df = - df/K;

end