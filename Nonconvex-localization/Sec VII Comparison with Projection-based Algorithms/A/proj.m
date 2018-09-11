function solution=proj(s,a,r)
% k=(s(2)-a(2))/(s(1)-a(1));
% b=(a(2)*s(1)-a(1)*s(2))/(s(1)-a(1));
% x = fsolve(@(x) (k*x+b-a(2))^2+(x-a(1))^2-r^2,s(1),optimset('Display','off'));
% y = k*x+b;
% solution=[x;y];
solution=a+r*(s-a)/norm(s-a);
end