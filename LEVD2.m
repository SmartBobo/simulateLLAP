function [output,newDC,Lmax,Lmin,freCount] = LEVD2(input,preDC,freCount)
%对输入的矩阵进行LEVD处理

Lmax = max(input);
Lmin = min(input);
s = zeros(120,1); 
temp = zeros(120,1);
%get variance
for b = 2:120
    temp(b) = input(b)-input(1); 
end
dsum = sum(temp)/120;
vsum = sum(sqrt(temp))/120;
if vsum + dsum * dsum > 15000  %15000 is power threshold
    if Lmax-Lmin > 220 * 4 %amplitude threshold
        newDC = ( 1 - 0.25)*preDC + 0.25*(Lmax+Lmin)/2;
    else
        newDC = 0;
    end
    freCount = 1;
else
    newDC = 0;
    
end
for t = 1:120
    s(t) = newDC;
end
output = s;
end