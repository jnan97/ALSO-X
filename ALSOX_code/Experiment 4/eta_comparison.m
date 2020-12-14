% Figure 6
clear;
clc;    
epsilon = 0.05;
upper_theta=0.05;
q=2;
k=0;
for theta =  0.00:0.01:upper_theta
k=k+1;

f=@(x,t) exp(-1/2*t.^2)/sqrt(2*pi).*(x-t).^q;
obj = @(x) integral(@(t) f(x,t),norminv(1-epsilon,0,1),x) - theta^q ;

U=theta/epsilon*1000+1000;    
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
result_2(k)=root-norminv(1-epsilon,0,1);

cvar_2(k)= normpdf(norminv(1-epsilon,0,1))/epsilon +theta/epsilon^(1/q)-norminv(1-epsilon,0,1);

end
q=3;
k=0;
for theta = 0.00:0.01:upper_theta
k=k+1;



f=@(x,t) exp(-1/2*t.^2)/sqrt(2*pi).*(x-t).^q;
obj = @(x) integral(@(t) f(x,t),norminv(1-epsilon,0,1),x) - theta^q ;

U=theta/epsilon*1000+1000;    
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
result_3(k)=root-norminv(1-epsilon,0,1);

cvar_3(k)= normpdf(norminv(1-epsilon,0,1))/epsilon +theta/epsilon^(1/q)-norminv(1-epsilon,0,1);

end


q=5;
k=0;
for theta = 0.00:0.01:upper_theta
k=k+1;


f=@(x,t) exp(-1/2*t.^2)/sqrt(2*pi).*(x-t).^q;
obj = @(x) integral(@(t) f(x,t),norminv(1-epsilon,0,1),x) - theta^q ;

U=theta/epsilon*1000+1000;    
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

result_5(k)=root-norminv(1-epsilon,0,1);

cvar_5(k)= normpdf(norminv(1-epsilon,0,1))/epsilon +theta/epsilon^(1/q)-norminv(1-epsilon,0,1);

end


q=10;
k=0;
for theta =  0.00:0.01:upper_theta
k=k+1;


f=@(x,t) exp(-1/2*t.^2)/sqrt(2*pi).*(x-t).^q;
obj = @(x) integral(@(t) f(x,t),norminv(1-epsilon,0,1),x) - theta^q ;

U=theta/epsilon*1000+1000;    
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
result_10(k)=root-norminv(1-epsilon,0,1);

cvar_10(k)= normpdf(norminv(1-epsilon,0,1))/epsilon +theta/epsilon^(1/q)-norminv(1-epsilon,0,1);

end

difference_2 = cvar_2-result_2;
difference_3 = cvar_3-result_3;
difference_5 = cvar_5-result_5;
difference_10 = cvar_10-result_10;

xx = linspace(0.00,upper_theta,upper_theta/0.01+1);


a=[difference_2;difference_3 ;difference_5;difference_10]';

xx_label = {'q=2','q=3','q=5','q=10'};

 plot(a(1,:), 'b-*','LineWidth',2,'MarkerSize',10)
 hold on 
 plot(a(2,:), 'r-x','LineWidth',2,'MarkerSize',10)
 hold on 
 plot(a(3,:), '-o','Color',[.61 .51 .74],'LineWidth',2,'MarkerSize',10)
 hold on 
 plot(a(4,:), '-pentagram','Color',[.71 .31 .94],'LineWidth',2,'MarkerSize',10)
 hold on 
  plot(a(5,:), 'm-hexagram','LineWidth',2,'MarkerSize',10)
 hold on 
  plot(a(6,:), '-diamond','Color',[.51 .21 .34],'LineWidth',2,'MarkerSize',10)
 hold on 
 
 yticks([0.20 0.30 0.40 0.50 0.60])
 set(gca,'XTick',1:4,'XTickLabel',xx_label,'FontSize',12)
 
 str1=' $$\theta=0.00$$';
  str2=' $$\theta=0.01$$';
   str3=' $$\theta=0.02$$';
    str4=' $$\theta=0.03$$';
    str5=' $$\theta=0.04$$';
    str6=' $$\theta=0.05$$';
 
 legend({str1,str2,str3,str4,str5,str6},'interpreter','latex','fontsize',16,'Location','northwest')
 xlabel('q-Wasserstein Ambiguity Set','Fontsize',18)
 ylabel('${\eta}^{\mathrm{CVaR}(q)}-{\eta}^{\mathrm{DRCCP}(q)}$','interpreter','latex','Fontsize',18)
% ylabel(' {\eta}^{\mathrm{DRCCP}(q)}','Fontsize',18);
 % title('10-Wasserstein','Fontsize',18)
 ylim([0.15 0.7])


