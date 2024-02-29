clc
clear
close all
gamma_bardB = 0:30;
gamma_bar = 10.^(gamma_bardB./10);
gammabar_SR = gamma_bar;
gammabar_RD = gamma_bar;
gammabar_hIdb = 1;
gammabar_hI = 10.^(gammabar_hIdb./10);
m_SR = 1;
L = 2;
N=10;
gamma_0db = 3;
gamma_0 = 10.^(gamma_0db./10);
s = (N *pi) / 4;
SIGMA = sqrt(N * (1 - (pi^2 / 16)));
C = ((s^2)/(2*SIGMA^2));
b_SR = 0.063;
omega_SR = 0.0007;
alpha_SR = ((2 * b_SR * m_SR) / ((2 * b_SR * m_SR) + omega_SR)^(m_SR))/ (2 * b_SR);

Delta_SR=[];
SOP = zeros(size(gamma_bar));
H = zeros(size(gamma_bar));

for index = 1:length(gamma_bar)
    delta_SR = omega_SR / (2 * b_SR * (2 * b_SR * m_SR + omega_SR));
            beta_SR = 1 / (2 * b_SR);
             Delta_SR(index) = (beta_SR - delta_SR) / gammabar_RD(index);
             B = (-s/(2*SIGMA^2*gammabar_RD(index)));
%             A = (1/(2*SIGMA^2*gammabar_RD(index)^2));
             A = 1;

    J_1 = [];
    for k_1 = 0:(m_SR-1)
        for t=0:k_1
            I_2=[];
            for v = 0:t
                I_1=[];
                for q = 0:(t-v)
                    I_1(q+1) = (nchoosek((t-v),q)*(gamma_0^(t-q))*(factorial(t+L-v-1)))/((exp(Delta_SR(index)*gamma_0))*(((1/gammabar_hI)+(Delta_SR(index)*(1+gamma_0)))^(t+L-v)));
                end
                I_2(v+1) = nchoosek(t,v)*sum(I_1);
            end
            J_1(k_1+1,t+1) = ((alpha_SR*pochhammer((1-m_SR),k_1) * ((-delta_SR)^k_1))/(((factorial(k_1))*factorial(t)*(gammabar_SR(index).^(k_1-t+1))*(Delta_SR(index).^(k_1-t+1))*(gammabar_hI^(L))*(factorial(L-1)))))*sum(I_2);
        end
    end
    
    J_2 = [];
    for k_1 = 0:(m_SR-1)
        for t = 0:k_1
            for k = 1:t
                Z_1=[];
                for p = 0:(2*k-1)
                    Z_2=[];
                    for v = 0:(t+p)
                        Z_3=[];
                        for q = 0:(t+p-v)
                            Z_3(q+1) = (nchoosek((t+p-v),q)*(gamma_0^(t+p-q))*(factorial(t+p+L-v)))/ (exp(Delta_SR(index)*gamma_0)*((1/gammabar_hI)+(Delta_SR(index)*(1+gamma_0)))^(t+p+L-v));
                        end
                        Z_2(v+1) = nchoosek((t+p),v)*sum(Z_3);
                    end
                    Z_1(p+1) = ((nchoosek((2*k-1),p)*B^(2*k-1-p)*(-1)^(k+1))/ (A^(k-0.5-p)*(2*k-1)*factorial(k-1)))*sum(Z_2);
                end
                J_2(k_1+1, t+1, k+1) = ((alpha_SR*pochhammer((1-m_SR),k_1) *((-delta_SR)^k_1)*(exp(((B^2)-(A*C))/A))) / ((sqrt(2*pi*A*(SIGMA^2)*(gammabar_RD(index)^2)))*(factorial(k_1))*(factorial(t))*(gammabar_SR(index)^(k_1+1))*(Delta_SR(index)^(k_1-t+1))*(gammabar_hI^(L))*(factorial(L-1))))*sum(Z_1);
            end
        end
    end
                        
    H1=[];
    for k = 1:t
        H1(k+1)=((((-1)^(k+1))^(B^(2*k-1)))/((A^(k-0.5))*(factorial(k-1))*(2*k-1)));
    end
    H(index) = ((1/(sqrt(2*pi*A*(SIGMA^2)*gammabar_RD(index)^2)))*(exp(((B^2)-(A*C))/A)))*sum(H1);

    SOP(index) = 1- (sum(sum(J_1))) - (H(index)*(sum(sum(J_1)))) + (sum(sum(sum(J_2))));

end
semilogy(gamma_bardB, SOP)  
xlabel('gamma-bar')
ylabel('SOP')
