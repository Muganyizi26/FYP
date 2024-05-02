clc
clear all
close all

gamma_bardB = 10:30;
gamma_bar = 10.^(gamma_bardB./10);
Lambdabar_Q = gamma_bar;
gammabar_SR = gamma_bar;
gammabar_RD = gamma_bar;
m_SR = 5;
m_Q = m_SR;
N=30;
gamma_0db = 2;
gamma_0 = 10.^(gamma_0db./10);
s = (N *pi) / 4;
SIGMA = sqrt(N * (1 - (pi^2 / 16)));
C = ((s^2)/(2*SIGMA^2));
b_SR = 0.251;
omega_SR = 0.279;
alpha_SR = ((2 * b_SR * m_SR) / ((2 * b_SR * m_SR) + omega_SR)^(m_SR))/ (2 * b_SR);
alpha_Q = alpha_SR;

Delta_SR=[]; Delta_Q=[];
B=[];
SOP = zeros(size(gamma_bar));
H = zeros(size(gamma_bar));
for index = 1:length(gamma_bar)
    delta_SR = omega_SR / (2 * b_SR * (2 * b_SR * m_SR + omega_SR));
    delta_Q = delta_SR;
            beta_SR = 1 / (2 * b_SR);
             Delta_SR(index) = (beta_SR - delta_SR) / gammabar_RD(index);
             Delta_Q = Delta_SR;
             B(index) = (-s/(2*SIGMA^2*gammabar_RD(index)));
             A = 1;
    J_1 =[];
    for k_1 = 0:(m_SR-1)
        for t = 0:k_1
            I_1 =[];
            for k_2 = 0:(m_Q-1)
                
                I_2=[];
                for v = 0:t
                    I_3 = [];
                    for q = 0:(t-v)
                        I_3(q+1) = (((nchoosek((t-v),q))*(gamma_0^(t-q))*factorial(t-v+k_2))/((exp(Delta_SR(index)*gamma_0))*(Delta_Q(index)+(Delta_SR(index)*(1+gamma_0)))^(t-v+k_2+1)));
                    end
                    I_2(v+1) = nchoosek(t,v)*sum(I_3);
                end
                I_1(k_2+1) = (((-delta_Q)^k_2)*pochhammer((1-m_Q),k_2))/((factorial(k_2))*(Lambdabar_Q(index))^(k_2+1))*sum(I_2);
            end
            J_1(k_1+1, t+1) = ((alpha_SR*pochhammer((1-m_SR),k_1) *((-delta_SR)^k_1))/((factorial(k_1))*(factorial(t))*(gammabar_SR(index)^(k_1+1))*(Delta_SR(index)^(k_1-t+1))))*alpha_Q*sum(I_1);
        end
    end
    
    J_2 = [];
    for k_1 = 0:(m_SR-1)
        for t =0:k_1
            for k = 1:t
                for p = 0:(2*k-1)
                    Z_1 = [];
                    for k_2 = 0:(m_Q-1)
                        Z_2 = [];
                        for v = 0:(t+p)
                            Z_3 = [];
                            for q = 0:(t+p-v)
                                Z_3(q+1) = ((nchoosek((t+p-v),q)*(gamma_0^(t+p-v))*(factorial(t+p-v+k_2)))/((exp(Delta_SR(index)*gamma_0))*(Delta_Q(index)+Delta_SR(index)*(1+gamma_0))^(t+p-v+k_2+1)));
                            end
                            Z_2(v+1) = nchoosek((t+p),v)*sum(Z_3);
                        end
                        Z_1(k_2+1) = ((((-delta_Q)^k_2)*pochhammer((1-m_Q),k_2))/((factorial(k_2))*(Lambdabar_Q(index)^(k_2+1))))*sum(Z_2);
                    end
                    J_2(k_1+1, t+1,k+1,p+1) = (((nchoosek((2*k-1),p))*(B(index)^(2*k-1-p))*((-1)^(k+1))*alpha_SR*pochhammer((1-m_SR),k_1)*((-delta_SR)^k_1)*(exp(((B(index)^2)-(A*C))/A)))...
                        /((A^(k-0.5-p))*(2*k-1)*(factorial(k-1))*(sqrt(2*pi*A*(SIGMA^2)*(gammabar_RD(index)^2)))*(factorial(k_1))*(factorial(t))*(gammabar_SR(index)^(k_1+1))*(delta_SR^(k_1-t+1))))...
                        *alpha_Q*sum(Z_1);
                end
            end
        end
    end
    H1=[];
    for k = 1:t

        H1(k+1)=((((-1)^(k+1))^(B(index)^(2*k-1)))/((A^(k-0.5))*(factorial(k-1))*(2*k-1)));
    end
    H(index) = ((1/(sqrt(2*pi*A*(SIGMA^2)*gammabar_RD(index)^2)))*(exp(((B(index)^2)-(A*C))/A)))*sum(H1);
    
    SOP(index) = 1- (sum(sum(J_1))) - (H(index)*(sum(sum(J_1)))) + (sum(sum(sum(sum(J_2)))));
end
semilogy(gamma_bardB, SOP)  
xlabel('γ̄ (dB)')
ylabel('Pr(Cs ≤ Co)')
title('SOP vs γ̄ ')