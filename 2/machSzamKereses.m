function reziduum = machSzamKereses(M1,kappa,pt,p0,lambda,L,D)
p1TartalyCso = (1 + (kappa-1)/2*M1^2)^(kappa/(1-kappa))*pt;

M2 = machSzamCsoVegen(M1,kappa,lambda,L,D);

p1pCsillag = 1/M1 * sqrt((kappa+1)/(2 + (kappa-1)*M1^2));
p2pCsillag = 1/M2 * sqrt((kappa+1)/(2 + (kappa-1)*M2^2));
p1Cso = p1pCsillag/p2pCsillag*p0;

reziduum = (p1Cso-p1TartalyCso)/p1TartalyCso;

end