function [z] = errFuncSparse(approx,f)
%Calculating \| approx - f\|_{L^2}.

%%% THIS CODE IS SLOW ------
%%%   -- We compute quadrature ponts on each element => slow
n = size(approx,1);

degree = 10;

z = 0;
for i=1:n
    for j=1:n
        a = 0+(j-1)*(1/n);
        b = a+1/n;
        d = 1-(i-1)*(1/n);
        c = d-(1/n);
        q = approx(i,j);
        
        
        [x,wx] = lgwt(degree,a,b);
        [y,wy] = lgwt(degree,c,d);
        
        [XX,YY] = meshgrid(x,y);
        
        %z = z + (1/9)*(a^3-b^3)*(c^3-d^3) ... %int_T (x*y-q)^2
        %      - (1/2)*q*(a^2-b^2)*(c^2-d^2) ...
        %      + q^2*(a-b)*(c-d);
        intf = sum( (f(XX,YY)-q).^2 .* (wx*wy') ,'all');
        z = z + intf;
    end
end
z = sqrt(z);

end