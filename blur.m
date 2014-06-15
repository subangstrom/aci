
function b=blur(a,w,window)
p=zeros(window,1);
for i=1:1:length(p)
    p(i)=exp(-(i-(window+1)/2)^2/w^2);
end
b=conv(a,p,'same');
end