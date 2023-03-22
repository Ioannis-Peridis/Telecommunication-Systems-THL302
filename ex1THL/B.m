%B.1 
%intialization of parameter A and vector of integer values k
A = 5;
k = 0:2*A;
%creation of srrc pulses
[f1,t1]=srrc_pulse(T,over,A,a(1));
[f2,t2]=srrc_pulse(T,over,A,a(2));
[f3,t3]=srrc_pulse(T,over,A,a(3));
%initialization of integrals to vectos filled with zeros
integral1=zeros(1,length(k));
integral2=integral1;
integral3=integral1;
%for loop for every value of vector k
for ki=0:length(k)-1
   %creation of the delayed srrc pulse for a=0
    f1kt = [zeros(1, over*ki) f1(1:end-ki*over)];
   %creation of product (srrc pulse)*(srrc pulse delayed by kT) for a=0
    product1 = f1.*f1kt;
    if(ki<4)
        disp('press any key to see the graphs of (srrc pulse) and (srrc pulse delayed by kT)for a=0:')
        pause
        figure()
        plot(t1,f1)
        hold on
        plot(t1,f1kt);
        title('for a=0:(srrc pulse) and (srrc pulse delayed by kT)')
        hold off
        disp('press any key to see the graph of(srrc pulse)*(srrc pulse delayed by kT) for a=0:')
        pause
        figure()
        plot(t1,product1)
        title('for a=0:(srrc pulse)*(srrc pulse delayed by kT)')
        %creation of integral of (srrc pulse)*(srrc pulse delayed by kT)
        %for a=0
        integral1(ki+1)=sum(product1)*Ts;
    end
   %creation of the delayed srrc pulse for a=0.5
   f2kt = [zeros(1, over*ki) f2(1:end-ki*over)];
   %creation of product (srrc pulse)*(srrc pulse delayed by kT) for a=0.5
    product2 = f2.*f2kt;
    if(ki<4)
        disp('press any key to see the graphs of (srrc pulse) and (srrc pulse delayed by kT)for a=0.5:')
        pause
        figure()
        plot(t2,f2)
        hold on
        plot(t2,f2kt);
        title('for a=0.5:(srrc pulse) and (srrc pulse delayed by kT)')
        hold off
        disp('press any key to see the graph of(srrc pulse)*(srrc pulse delayed by kT) for a=0.5:')
        pause
        figure()
        plot(t2,product2)
        title('for a=0.5:(srrc pulse)*(srrc pulse delayed by kT)')
        %creation of integral of (srrc pulse)*(srrc pulse delayed by kT)
        %for a=0.5
        integral2(ki+1)=sum(product2)*Ts;
    end
    %creation of the delayed srrc pulse for a=1
    f3kt = [zeros(1, over*ki) f3(1:end-ki*over)];
    %creation of product (srrc pulse)*(srrc pulse delayed by kT) for a=1
    product3 = f3.*f3kt;
    if(ki<4)
        disp('press any key to see the graphs of (srrc pulse) and (srrc pulse delayed by kT)for a=1:')
        pause
        figure()
        plot(t3,f3)
        hold on
        plot(t3,f3kt);
        title('for a=1:(srrc pulse) and (srrc pulse delayed by kT)')
        hold off
        disp('press any key to see the graph of(srrc pulse)*(srrc pulse delayed by kT) for a=1:')
        pause
        figure()
        plot(t3,product3)
        title('for a=1:(srrc pulse)*(srrc pulse delayed by kT)')
        %creation of integral of (srrc pulse)*(srrc pulse delayed by kT)
        %for a=1
        integral3(ki+1)=sum(product3)*Ts;
    end
end
disp(' ')
disp('the values of the integral of (srrc pulse)*(srrc pulse delayed by kT) for a=0,0.5,1 and k=0,1,2,3 ')
disp(' ')
integral1(1:4)
integral2(1:4)
integral3(1:4)