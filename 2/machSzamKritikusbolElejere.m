function reziduum = machSzamKritikusbolElejere(M1,kappa,lambda,xmax,D)
reziduum = (1 - M1^2)/(kappa*M1^2) + (kappa+1)/(2*kappa)*log(((kappa+1)*M1^2)/(2 + (kappa-1)*M1^2)) - lambda*xmax/D;
end