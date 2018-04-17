function [M2, xmax] = machSzamCsoVegen(M1,kappa,lambda,L,D)
lambdaXmaxD = (1 - M1^2)/(kappa*M1^2) + (kappa+1)/(2*kappa)*log(((kappa+1)*M1^2)/(2 + (kappa-1)*M1^2));
xmax = lambdaXmaxD/lambda*D;

L2 = xmax - L;

options = optimset('Display','off');
M2 = fsolve(@(M2)machSzamKritikusbolElejere(M2,kappa,lambda,L2,D),M1,options);
end
