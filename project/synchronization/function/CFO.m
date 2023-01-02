function result = CFO(x,fs,dw,phi)
%e^j(wt+phi)
N=0:length(x)-1;
t=N'/fs;
result=x.*exp(1i.*(dw.*t+phi));
end

