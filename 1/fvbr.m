function dy = fvbr(t,y,m,k,s,f0,w)

dy = zeros(2,1);
dy(1) = y(2);
dy(2) = 1/m*(-k*y(2)-s*y(1)+f0*sin(w*t));

end
