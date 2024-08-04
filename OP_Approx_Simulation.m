function [ff,xx]=P_out_sim_anal(N,eta,SNRdB,minvalue,maxvalue)
tic
n = 1e6;
%gamthdb =0:1:30;
SNR=10^(SNRdB/10);
Xm=zeros(1,n);
v = 0.2;
s = sqrt((1-v*v)/2);
% cpt=zeros(1,length(gamthdb));
 f = abs(randn(1,n)+1i*randn(1,n));
   
    for k=1:N
    hk = normrnd( 0, s , [1,n])+1i*normrnd( v, s , [1,n]);
    gk = normrnd( 0, s , [1,n])+1i*normrnd( v, s , [1,n]);
    h = abs(hk);
    g = abs(gk);
    Xk = h.*g;
    Xm = Xm + Xk;
    end

   Ym=Xm.*f;
    gam=10*log10(eta*SNR*Ym.^2);
    m1=mean(Xm)
    m2=mean(Xm.^2)
   m3=mean(Xm.^3)
    m4=mean(Xm.^4)
    [mu_1,mu_2,mu_3,mu_4] = moments(s,N,v);
psi2 = mu_2/mu_1;
psi3 = mu_3/mu_2;
psi4 = mu_4/mu_3;

a3 = (4*psi4 - 9*psi3 + 6*psi2 - mu_1)/(-psi4 + 3*psi3 - 3*psi2 + mu_1)
a2 = (a3/2)*(psi4 - 2*psi3 + psi2) + 2*psi4 -3*psi3 + psi2
a6 = ((a3*(psi2-mu_1)+2*psi2 - mu_1)/a2) - 3
a7 = sqrt(((((a3*(psi2-mu_1)+2*psi2 - mu_1)/a2) - 1)^2)-((4/a2)*mu_1*(a3+1)))
a4 = (a6+a7)/2
a5 = (a6-a7)/2
a1 = gammaz(a3+1)/(a2*gammaz(a4+1)*gammaz(a5+1))
    
%xRange = 0:30.1;                %# Range of integers to compute a probability for
% A = hist(Ym,xRange);        %# Bin the data
% plot(xRange,A./numel(Ym),'ro','LineWidth',1);  %# Plot the probabilities for each integer
% h=cdfplot(gam);
Lim=maxvalue-minvalue;
%esp=10^(-4);
[f, x] = ecdf(gam);
j=min(find(x>=minvalue))
i=1;
xx(i)=x(j);
ff(i)=f(j);
j=j+1;
i=i+1;
while (i<=Lim-1)
   while  x(j)<xx(i-1)+1
       j=j+1;
   end 
   xx(i)=x(j);
   ff(i)=f(j);
   i=i+1;
end
%  for i=1:Lim
%      if i<=20
%          j=1+(i-1)*inn/3
%         % j=1;
%    %  else
%     % j=1+(i-1)*inn/10;
%      %else
%      else
%      j=1+i*inn; 
%      end
%  
%  xx(i)=x(j)
%  ff(i)=f(j)
%  end
% e=0;
% i=1;
% while (i<=length(x)|| e~=1)
%     floor(x(i)*100
% end
% for i=1:Lim
% ind(i)=min(find(x>=i-esp/2 & x<i+eps/2));
% end
% xx=x(ind);
% ff=f(ind);
%plot(xx,ff,'r*')
%xx(1:2:end)
%spacebtwmarker=400000;
%xx=xx+30; %dBm
plot(xx,ff,'b*');
%plot(xx(1:spacebtwmarker:end),ff(1:spacebtwmarker:end),'r*')
% set( h, 'Color', 'r');
 xlabel('\gamma_{th}');
 ylabel('P_{out}');
%outageProb_sum_DREH(N,eta,SNRdB);
print -depsc fig1

timeElapsed = toc
end