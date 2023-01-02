function [corrected, error] = GardenerK(samples,k,M)
error = zeros(1,ceil(length(samples)/M));
corrected = zeros(size(error));
corrected(1) = samples(1);
for n=1:length(error)-1 %% n downsample symbol length
    interpolate = interp1(1:M+1,samples(M*(n-1)+1:M*n+1),[M/2+1 M+1]+M*error(n),'pchip');
    corrected(n+1) = interpolate(2); 
    %%   Corrected(1)= samples(1)  ; MidY=interpolate(1)=samples(M/2 +1) ; Corrected(2)= samples(M+1)
    error(n+1) = error(n) - 2*k*real(interpolate(1)*(conj(corrected(n+1)) - conj(corrected(n))));
end
corrected=corrected';
end


