function Q = qij_3state(i,j,n,lambda,delta) 
%This function computes the inverse laplace transforms in the case m=0; exact solutions are provided for all n \in \mathbb{N}. 
sum1=0; 
sum2=0;
sum3=0;
sum4=0; 
    
if abs(sum(lambda(1)+lambda(3))-sum(lambda(2)+lambda(4)))>1e-3 
    if i==0 && j==0
        for i=1:(n+1) 
        a = ((-1)^(n-i+1))*nchoosek_ln(2*n-i,n-i+1)/((lambda(2)+lambda(4)-lambda(3)-lambda(1))^(2*n+1-i)); 
        sum2 = sum2 + a*(delta^(i-1))*exp(-(lambda(1)+lambda(3))*delta)/factorial(i-1);
        end 

        for i=1:n 
        b = ((-1)^(n-i))*nchoosek_ln(2*n-i,n-i)/((lambda(1)+lambda(3)-lambda(2)-lambda(4))^(2*n+1-i)); 
        sum1 = sum1 + b*(delta^(i-1))*exp(-(lambda(2)+lambda(4))*delta)/factorial(i-1); 
        end  

        Q = ((lambda(1)*lambda(2))^n)*(sum1+sum2); 

    elseif i==0 && j==1 
        for i=1:n 
        a = ((-1)^(n-i))*nchoosek_ln(2*n-i-1,n-i)/((lambda(2)+lambda(4)-lambda(3)-lambda(1))^(2*n-i)); 
        sum2 = sum2 + a*(delta^(i-1))*exp(-(lambda(1)+lambda(3))*delta)/factorial(i-1);
        end 

        for i=1:n 
        b = ((-1)^(n-i))*nchoosek_ln(2*n-i-1,n-i)/((lambda(1)+lambda(3)-lambda(2)-lambda(4))^(2*n-i)); 
        sum1 = sum1 + b*(delta^(i-1))*exp(-(lambda(2)+lambda(4))*delta)/factorial(i-1); 
        end 

        Q = ((lambda(1))^n)*((lambda(2))^(n-1))*(sum1+sum2); 

    elseif i==1 && j==1 
        for i=1:n 
        a = ((-1)^(n-i))*nchoosek_ln(2*n-i,n-i)/((lambda(2)+lambda(4)-lambda(3)-lambda(1))^(2*n-i+1)); 
        sum2 = sum2 + a*(delta^(i-1))*exp(-(lambda(1)+lambda(3))*delta)/factorial(i-1);
        end 

        for i=1:(n+1) 
        b = ((-1)^(n-i+1))*nchoosek_ln(2*n-i,n-i+1)/((lambda(1)+lambda(3)-lambda(2)-lambda(4))^(2*n-i+1)); 
        sum1 = sum1 + b*(delta^(i-1))*exp(-(lambda(2)+lambda(4))*delta)/factorial(i-1);  
        end 

        Q = ((lambda(1)*lambda(2))^n)*(sum1+sum2); 

    elseif i==1 && j==0 
        for i=1:(n+1) 
        a = ((-1)^(n-i+1))*nchoosek_ln(2*n-i+1,n-i+1)/((lambda(2)+lambda(4)-lambda(3)-lambda(1))^(2*n-i+2)); 
        sum2 = sum2 + a*(delta^(i-1))*exp(-(lambda(1)+lambda(3))*delta)/factorial(i-1);
        end 

        for i=1:(n+1) 
        b = ((-1)^(n-i+1))*nchoosek_ln(2*n-i+1,n-i+1)/((lambda(1)+lambda(3)-lambda(2)-lambda(4))^(2*n-i+2)); 
        sum1 = sum1 + b*(delta^(i-1))*exp(-(lambda(2)+lambda(4))*delta)/factorial(i-1);  
        end 

        Q = ((lambda(1))^n)*((lambda(2))^(n+1))*(sum1+sum2);

    elseif i==0 && j==2     
        for i=1:(n+1) 
        a = ((-1)^(n-i+1))*nchoosek_ln(2*n-i,n-i+1)/((lambda(2)+lambda(4)-lambda(3)-lambda(1))^(2*n+1-i));
        p = poisscdf(i-1,(lambda(1)+lambda(3))*delta);
        sum1 = sum1 + a*(1-p)/(lambda(1)+lambda(3))^i;  
        end 

        for i=1:n 
        b = ((-1)^(n-i))*nchoosek_ln(2*n-i,n-i)/((lambda(1)+lambda(3)-lambda(2)-lambda(4))^(2*n+1-i)); 
        p = poisscdf(i-1,(lambda(2)+lambda(4))*delta);
        sum2 = sum2 + b*(1-p)/(lambda(2)+lambda(4))^i;  
        end  

        for i=1:n 
        a = ((-1)^(n-i))*nchoosek_ln(2*n-i-1,n-i)/((lambda(2)+lambda(4)-lambda(3)-lambda(1))^(2*n-i)); 
        p = poisscdf(i-1,(lambda(1)+lambda(3))*delta);
        sum3 = sum3 + a*(1-p)/(lambda(1)+lambda(3))^i;  
        end 

        for i=1:n 
        b = ((-1)^(n-i))*nchoosek_ln(2*n-i-1,n-i)/((lambda(1)+lambda(3)-lambda(2)-lambda(4))^(2*n-i)); 
        p = poisscdf(i-1,(lambda(2)+lambda(4))*delta);
        sum4 = sum4 + b*(1-p)/(lambda(2)+lambda(4))^i;
        end 

        Q = lambda(3)*((lambda(1)*lambda(2))^n)*(sum1+sum2) + lambda(4)*((lambda(1)*lambda(2))^n)*(sum3+sum4)/lambda(2);

    elseif i==1 && j==2 
        for i=1:(n+1) 
        a = ((-1)^(n-i+1))*nchoosek_ln(2*n-i+1,n-i+1)/((lambda(2)+lambda(4)-lambda(3)-lambda(1))^(2*n-i+2)); 
        p = poisscdf(i-1,(lambda(1)+lambda(3))*delta);
        sum1 = sum1 + a*(1-p)/(lambda(1)+lambda(3))^i;  
        end 

        for i=1:(n+1)
        b = ((-1)^(n-i+1))*nchoosek_ln(2*n-i+1,n-i+1)/((lambda(1)+lambda(3)-lambda(2)-lambda(4))^(2*n-i+2)); 
        p = poisscdf(i-1,(lambda(2)+lambda(4))*delta);
        sum2 = sum2 + b*(1-p)/(lambda(2)+lambda(4))^i;  
        end  

        for i=1:n 
        a = ((-1)^(n-i))*nchoosek_ln(2*n-i,n-i)/((lambda(2)+lambda(4)-lambda(3)-lambda(1))^(2*n-i+1)); 
        p = poisscdf(i-1,(lambda(1)+lambda(3))*delta);
        sum3 = sum3 + a*(1-p)/(lambda(1)+lambda(3))^i;
        end 

        for i=1:(n+1) 
        b = ((-1)^(n-i+1))*nchoosek_ln(2*n-i,n-i+1)/((lambda(1)+lambda(3)-lambda(2)-lambda(4))^(2*n-i+1)); 
        p = poisscdf(i-1,(lambda(2)+lambda(4))*delta);
        sum4 = sum4 + b*(1-p)/(lambda(2)+lambda(4))^i;
        end 

        Q = lambda(3)*(lambda(1)^n)*(lambda(2)^(n+1))*(sum1+sum2) + lambda(4)*((lambda(1)*lambda(2))^n)*(sum3+sum4);  
    end 
else 
    if i==0 && j==0 
        if n==0 
            Q = exp(-(lambda(2)+lambda(4))*delta);
        else 
            Q = lambda(1)^n * lambda(2)^(n+i-j) * delta^(2*n+i-j) * exp(-(lambda(2)+lambda(4))*delta)/factorial(2*n+i-j);
        end 
    elseif i==0 && j==1 
        if n==0 
            Q = 0; 
        else 
            Q = lambda(1)^n * lambda(2)^(n+i-j) * delta^(2*n+i-j) * exp(-(lambda(2)+lambda(4))*delta)/factorial(2*n+i-j);
        end 
    elseif i==1 && j==0 
        if n==0 
            Q = delta*lambda(2)*exp(-(lambda(2)+lambda(4))*delta);
        else 
            Q = lambda(1)^n * lambda(2)^(n+i-j) * delta^(2*n+i-j) * exp(-(lambda(2)+lambda(4))*delta)/factorial(2*n+i-j);
        end 
    elseif i==1 && i==1 
        if n==0 
            Q = exp(-(lambda(2)+lambda(4))*delta);
        else 
            Q = lambda(1)^n * lambda(2)^(n+i-j) * delta^(2*n+i-j) * exp(-(lambda(2)+lambda(4))*delta)/factorial(2*n+i-j);
        end 
    elseif i==0 && j==2 
        if n==0 
            Q = lambda(3)*(1- exp(-(lambda(2)+lambda(4))*delta))/lambda(2)+lambda(4);
        else 
            p1 = poisscdf(2*n,(lambda(2)+lambda(4))*delta);
            p2 = poisscdf(2*n-1,(lambda(2)+lambda(4))*delta);
            Q = (lambda(3)*lambda(1)^n * lambda(2)^n*(1-p1)/(lambda(2)+lambda(4))^(2*n+1)) + (lambda(4)*lambda(1)^n * lambda(2)^(n-1)*(1-p2)/(lambda(2)+lambda(4))^(2*n));
        end 
    elseif i==1 && j==2 
        if n==0 
            Q = (lambda(3)*lambda(2)*(1-exp(-(lambda(2)+lambda(4))*delta)*(1-(lambda(2)+lambda(4))*delta))/(lambda(2)+lambda(4))^2) +lambda(4)*(1- exp(-(lambda(2)+lambda(4))*delta))/lambda(2)+lambda(4);
        else 
            p1 = poisscdf(2*n+1,(lambda(2)+lambda(4))*delta);
            p2 = poisscdf(2*n,(lambda(2)+lambda(4))*delta);
            Q = (lambda(3)*lambda(1)^n * lambda(2)^(n+1)*(1-p1)/(lambda(2)+lambda(4))^(2*n+2)) + (lambda(4)*lambda(1)^n * lambda(2)^(n)*(1-p2)/(lambda(2)+lambda(4))^(2*n+1));
        end 
    end 
end 



    
