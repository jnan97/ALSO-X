% Figure 5
clear;
clc;    
theta = 0.20;
k=0;
q=2;
for epsilon = 0.01:0.01:0.10
k=k+1;

f=@(x,t) exp(-1/2*t.^2)/sqrt(2*pi).*(x-t).^q;
obj = @(x) integral(@(t) f(x,t),norminv(1-epsilon,0,1),x) - theta^q ;


U=theta/epsilon*1000;    
L=norminv(1-epsilon,0,1);

while U-L>1e-10    
    root=(U+L)/2;    
    if obj(root)==0    
        break;    
    end
    if obj(root)*obj(U)<0    
        L=root;
    else
        U=root;
    end
end
eta_q(k)=root;

drccp_q(k)=root-norminv(1-epsilon,0,1);

cvar(k)= normpdf(norminv(1-epsilon,0,1))/epsilon +theta/epsilon^(1/q)-norminv(1-epsilon,0,1);
end

k=0;
for epsilon = 0.01:0.01:0.10
k=k+1;

g = @(alpha) alpha +((1-normcdf(alpha,0,1))^(-1/q)*theta*(q-1)/q) - eta_q(k);
U=eta_q(k);    
L=0;
while U-L>1e-10    
    root=(U+L)/2;    
    if g(root)==0    
        break;    
    end
    if g(root)*g(U)<=0    
        L=root;
    else
        U=root;
    end
end
alpha_hat(k)=root;
end

k=0;
for epsilon = 0.01:0.01:0.10
k=k+1;
f=@(x) theta/q*(1-normcdf(alpha_hat(k),0,1))^((q-1)/q) + normpdf(alpha_hat(k),0,1)-alpha_hat(k)+alpha_hat(k)* normcdf(alpha_hat(k),0,1)-(normpdf(eta_q(k)-x,0,1)-(eta_q(k)-x)+(eta_q(k)-x)* normcdf(eta_q(k)-x,0,1));

U=theta/epsilon*1000+1000;     
L=0;

while U-L>1e-10    
    root=(U+L)/2;    
    if f(root)==0    
        break;    
    end
    if f(root)*f(U)<0    
        L=root;
    else
        U=root;
    end
end
result_alsox(k)=root;
end

xx = linspace(0.01,0.10,10);

yy = theta ./ xx.^(1/q);

plot(xx,result_alsox,'r-s','LineWidth',2,'MarkerSize',10)
hold on
plot(xx,cvar, 'm-*','LineWidth',2,'MarkerSize',10)
hold on 
plot(xx,drccp_q,'b->','LineWidth',2,'MarkerSize',10)
hold on 
plot(xx,yy,'k-x','LineWidth',2,'MarkerSize',10)
str1=' $${\eta}^{\mathrm{ALSO-X}(2)}$$';
str2=' $${\eta}^{\mathrm{CVaR}(2)}$$';
str3=' $${\eta}^{\mathrm{DRCCP}(2)}$$';
str4=' $${\eta}^{\mathrm{VaR}(2)}$$';
legend({str1,str2,str3,str4},'interpreter','latex','fontsize',18,'Location','northeast')
xlabel('Risk Level \epsilon','Fontsize',18)
ylabel('Level of Robustness \eta','Fontsize',18);
ylim([0.2 2.8])
xlim([0.01 0.10])
%  xlim([0.05 0.50])
%title('\eta comparisons')
%xlabel('different \theta when \epsilon = 0.05')