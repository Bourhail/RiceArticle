function xxx=PDF_sim_anal(N)
n = 10^6;
Xm = zeros(1,n);
v = 0.25;
s = sqrt((1-v*v)/2);
r2 = abs(normrnd( 0, s , [1,n])+1i*normrnd( v, s , [1,n]));
    for k=1:N
    hk = normrnd( 0, s , [1,n])+1i*normrnd( v, s , [1,n]);
    gk = normrnd( 0, s , [1,n])+1i*normrnd( v, s , [1,n]);
    h = abs(hk);
    g = abs(gk);
    Xk = h.*g;
    Xm = Xm + Xk;
    end
R = Xm.*r2;

mean(Xm);
mean(Xm.^2);
mean(Xm.^3);
mean(Xm.^4);
 xRange = 0:70;
 A = hist(R,xRange);
 y =0.1:0.1:70;
  W = 50;
 [mu_1,mu_2,mu_3,mu_4] = moments(s,N,v);
psi2 = mu_2/mu_1;
psi3 = mu_3/mu_2;
psi4 = mu_4/mu_3;

a3 = (4*psi4 - 9*psi3 + 6*psi2 - mu_1)/(-psi4 + 3*psi3 - 3*psi2 + mu_1);
a2 = (a3/2)*(psi4 - 2*psi3 + psi2) + 2*psi4 -3*psi3 + psi2;
a6 = ((a3*(psi2-mu_1)+2*psi2 - mu_1)/a2) - 3;
a7 = sqrt(((((a3*(psi2-mu_1)+2*psi2 - mu_1)/a2) - 1)^2)-((4/a2)*mu_1*(a3+1)));
a4 = (a6+a7)/2;
a5 = (a6-a7)/2;
a1 = gammaz(a3+1)/(a2*gammaz(a4+1)*gammaz(a5+1));
k_max = 1;
 for i=1:length(y)
     pdf_DR(i)=0;
     for k=0:k_max
   f=@(t) (a1./(sqrt(2).*s)).*(exp((-v*v)/(2*s*s))).*(1/(gammaz(k+1).*factorial(k))).*(((v*v)/(2*s*s))^(k)).*((gammaz(a4+t).*gammaz(a5+t).*gammaz((2*k+t+1)./2))./(gammaz(a3+t))).*((y(i)/(sqrt(2)*a2*s)).^(-t));
  cs=max(max(real(-a4),real(-a5)),real(-2*k-1))+1;
 b = real(1/(2*pi*1i)*integral(f,cs-W*1i,cs+W*1i));
 pdf_DR(i)=pdf_DR(i)+b;
%   pdf_DR(i) = a1*meijerG([],[a3],[a4,a5],[], y(i)/a2);
%   pdf_DR(1:length(y)) = a1*meijerG([],[a3],[a4,a5],[], y/a2);
     end
 end
 plot(xRange,A./numel(R),'bo',y,pdf_DR,'g-','LineWidth',1.5);
 xlabel('r');
ylabel('f_{R}(r)');
legend('Simulation','Approximation');
title('The PDF of the e2e cascaded channel gain')
end