clc
clear
close all
gamma_bardB = 0:30;
gamma_bar = 10.^(gamma_bardB./10);
gammabar_SR = gamma_bar;
% gammabar_SE = gamma_bar;
gammabar_RD = gamma_bar; %find out more even in the derivation for F_RD
Lambdabar_1 = gamma_bar;
gammabar_hIdb = 3; % SNR of the interference link
gammabar_hI = 10.^(gammabar_hIdb./10); 
m_SR = 1;
m_SE = m_SR;
L = 1; % Number of interferers 
N=30; % Number of RIS reflecting elements
gammabar_SEdB = 3;
gammabar_SE = 10.^(gammabar_SEdB./10); % SNR of the eavesdropper
gamma_0db = 2; % target secrecy threshold SNR
gamma_0 = 10^(gamma_0db/10);
s = (N *pi) / 4;
SIGMA = sqrt(N * (1 - (pi^2 / 16)));
C = ((s^2)/(2*SIGMA^2));
b_SR = 0.063;
omega_SR = 0.0007;
alpha_SR = ((2 * b_SR * m_SR) / ((2 * b_SR * m_SR) + omega_SR)^(m_SR))/ (2 * b_SR);
alpha_SE = alpha_SR;

Delta_SR=[];
Delta_SE=[];
B=[];
H = zeros(size(gamma_bar));
SOP = zeros(size(gamma_bar));
for index = 1:length(gamma_bar)
    delta_SR = omega_SR / (2 * b_SR * (2 * b_SR * m_SR + omega_SR));
    delta_SE = delta_SR;
            beta_SR = 1 / (2 * b_SR);
             Delta_SR(index) = (beta_SR - delta_SR) / gammabar_RD(index);
             Delta_SE(index) = Delta_SR(index);
             B(index) = (-s/(2*SIGMA^2.*gammabar_RD(index)));
%             A = (1/(2*SIGMA^2*gammabar_RD(index)^2));
             A = 1;
    G_1 = [];
    for k_2=0:(m_SR-1)
        I_1 = [];
        for p=0:k_2
            I_2=[];
            for k_1 = 0:(m_SE-1)
                I_3=[];
                for v=0:(k_1+1)
                    I_4=[];
                    for k=0:p
                        I_5=[];
                        for m=0:(p-k)
                            I_6=[];
                            for z=0:m
                                I_6(z+1)= (((-1)^z)*nchoosek((L+v+z-1),z)*(Delta_SE(index)^z)*(factorial(k_1+p-k+z)))/((Delta_SE(index) + Delta_SR(index)*(1+gamma_0))^(k_1+p-k+z-1));
                            end
                            I_5(m+1) = nchoosek((p-k),m)*(gamma_0^(p-m))*(gammabar_hI^(v+z))*sum(I_6);
                        end
                        I_4(k+1)= nchoosek(p,k)*sum(I_5);
                    end
                    I_3(v+1) = nchoosek((k_1+1),v)*(factorial(L+v-1))*exp(-Delta_SR(index)*gamma_0)*sum(I_4);
                end
                I_2(k_1+1) = ((((-delta_SE)^k_1)*pochhammer((1-m_SE),k_1))/(((factorial(k_1))^2)*(gammabar_SE^(k_1+1))*(factorial(L-1))))*sum(I_3);
            end
            I_1(p+1) = ((((-delta_SR)^k_2)*pochhammer((1-m_SR),k_2))/((factorial(k_2))*(gammabar_SR(index)^(k_2+1))*factorial(p)*(Delta_SR(index)^(k_2-p+1))))*alpha_SE*sum(I_2);
        end
        G_1(k_2+1) = alpha_SR*sum(I_1);
    end
    
    G_2 = [];
    for k_2 =0:(m_SR-1)
        E_1 =[];
        for p=0:k_2
            E_2 =[];
            for n=0:10
                E_3 =[];
                for r=0:(2*n+1)
                    E_4 =[];
                    for q=0:(2*n+1-r)
                        E_5 =[];
                        for k_1=0:(m_SE-1)
                            E_6 =[];
                            for v=0:(k_1+1)
                                E_7 =[];
                                for k=0:(p+r+q)
                                    E_8 = [];
                                    for m=0:(p+r+q-k)                                        
                                        E_9=[];
                                        for z=0:m
                                            E_9(z+1)= (((-1)^z)*nchoosek((L+v+z-1),z)*(Delta_SE(index)^z)*(factorial(k_1+p+r+q-k+z)))/((Delta_SE(index) + Delta_SR(index)*(1+gamma_0))^(k_1+p+r+q-k+z-1));
                                        end
                                        E_8(m+1) = nchoosek((p+r+q-k),m)*(gamma_0^(p+r+q-m))*(gammabar_hI^(v+z))*sum(E_9);
                                    end
                                    E_7(k+1)= nchoosek((p+r+q),k)*sum(E_8);
                                end
                                E_6(v+1) = nchoosek((k_1+1),v)*(factorial(L+v-1))*exp(-Delta_SR(index)*gamma_0)*sum(E_7);
                            end
                            E_5(k_1+1) = (((-delta_SE^(k_1))*pochhammer((1-m_SE),k_1))/(((factorial(k_1))^2)*(gammabar_SE^(k_1+1))*(factorial(L-1))))*sum(E_6);
                        end
                        E_4(q+1) = nchoosek((2*n+1-r),q)*((sqrt(A))^(r+q))*((B(index)/sqrt(A))^(2*n+1-r-q))*(factorial(L-1+r))*alpha_SE*sum(E_5);
                    end
                    E_3(r+1) = nchoosek((2*n+1),r)*sum(E_4);
                end
                E_2(n+1) = (((-1)^n)/((factorial(n))*(2*n+1)))*sum(E_3);
            end
            E_1(p+1) = ((((-delta_SR)^k_2)*pochhammer((1-m_SR),k_2))/((factorial(k_2))*(gammabar_SR(index)^(k_2+1))*(factorial(p))*(Delta_SR(index)^(k_2-p+1))*(sqrt(pi))*(factorial(L-1))))*sum(E_2);
        end
        G_2 = alpha_SR*sum(E_1);
    end
    
    H(index) = 0.5*erf(s/(sqrt(2*SIGMA^2*Lambdabar_1(index)^2)));
    SOP(index) = 1-(sum(G_1))-(H(index))*(sum(G_1))-(sum(G_2));

end
semilogy(gamma_bardB, SOP)
xlabel('Average SNR γ̄ (dB)')
ylabel('Secrecy Outage probability')
title('SOP vs γ̄ ')

