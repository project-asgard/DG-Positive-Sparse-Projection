function [proj] = calcProj(B,f)
%Calculate L2 projection onto sparse basis \hat{V}_N^0.

M = numel(B);

n = size(B{1},1);
proj = zeros(n);

degree = 10;

%Calculate integral of f over each fine mesh
F = zeros(n);
for i=1:n
    for j=1:n
        %calculate rectangle [a,b]x[c,d] coordinate
        a = 0+(j-1)*(1/n);
        b = a+1/n;
        d = 1-(i-1)*(1/n);
        c = d-(1/n);
        
        [x,wx] = lgwt(degree,a,b);
        [y,wy] = lgwt(degree,c,d);
        
        [XX,YY] = meshgrid(x,y);
        F(i,j) = sum(sum(f(XX,YY).*(wx*wy')));
    end
end

for m=1:M
    proj = proj + sum(sum(B{m}.*F))*B{m};
end


end

