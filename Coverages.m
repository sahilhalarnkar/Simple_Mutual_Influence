%% Figure 1A
% Run this section to generate Figure 1A
thetaR1=[];
thetaR2=[];
thetaH=[];
thetastar=[];
r1=[];
r2=[];
r3=[];

% Initial point for the coverages of R1, R2 and H respectively 
x1=[1 1 1];


% Iterating over different relative concentrations of CR2
for CR2=0:0.1:4
    % Rate Parameters = [k2,k-2,k4,k-4,k1,k-1,k3,k5,CR1,CR2,CH]
    k1=[1,0.1,1,0.1,2,0.1,1,1,1,CR2,1];

    x1=lsqnonlin(@(s) ss_LH3(s,k1),[x1(1);x1(2);x1(3)],[0 0 0],[1 1 1]);
    if(CR2==0)
        x2=[x1(1) x1(3)];
    end
    
    % Storing reaction rates and coverages in lists
    r1=[r1 k1(7)*x1(1)*x1(3)];
    r2=[r2 k1(8)*x1(2)*x1(3)];
    r3=[r3 k1(7)*x2(1)*x2(2)];
    thetaR1=[thetaR1 x1(1)];
    thetaR2=[thetaR2 x1(2)];
    thetaH=[thetaH x1(3)];    

end

figure;
ax1=subplot(2,1,1);
plot(0.1:0.1:4,enhance(2:end));
ylabel('EF')
ax2=subplot(2,1,2);
plot(0:0.1:4,thetaR1,0:0.1:4,thetaR2,0:0.1:4,thetaH);
ylabel('Coverage')

ax1.Position(2)=ax2.Position(2)+ax2.Position(4);
linkaxes([ax1 ax2],'x')

%% Figure S1A
% Run this section to generate Figure S1A
thetaR1=[];
thetaR2=[];
thetaH=[];
thetastar=[];
r1=[];
r2=[];
r3=[];

% Initial point for the coverages of R1, R2 and H respectively 
x1=[1 1];


% Iterating over different relative concentrations of CR2
for CR2=0:0.1:4
    % Rate Parameters = [k2,k-2,k4,k-4,k1,k-1,k3,k5,CR1,CR2,CH]
    k1=[1,0.1,1,0.1,2,0.1,1,1,1,CR2,1];
    
    x1=lsqnonlin(@(s) ss_ER(s,k1),[x1(1);x1(2)],[0 0] ,[1 1]);
    if(CR2==0)
        x2=[x1(1) x1(2)];
    end
    
    % Storing reaction rates and coverages in lists
    r1=[r1 k1(7)*x1(1)*x1(2)];
    r2=[r2 k1(8)*x1(2)*x1(2)];
    r3=[r3 k1(7)*x2(1)*x2(2)];
    thetaR1=[thetaR1 x1(1)];
    thetaH=[thetaH x1(2)];    

end

enhance=r1./r3;
select=r1./r2;

figure;
ax1=subplot(2,1,1);
plot(0.1:0.1:4,enhance(2:end));
ylabel('EF')
ax2=subplot(2,1,2);
plot(0:0.1:4,thetaR1,0:0.1:4,thetaH);
ylabel('Coverage')
ax1.Position(2)=ax2.Position(2)+ax2.Position(4);
linkaxes([ax1 ax2],'x')
%% Figure S2A
% Run this section to generate Figure S2A
thetaR1=[];
thetaR2=[];
thetaH=[];
thetastar=[];
r1=[];
r2=[];
r3=[];

% Initial point for the coverages of R1, R2 and H respectively 
x1=[1 1 1];


% Iterating over different relative concentrations of CR2
for CR2=0:0.1:4
    % Rate Parameters = [k2,k-2,k4,k-4,k1,k-1,k3,k5,CR1,CR2,CH]
    k1=[1,0.1,1,0.1,2,0.1,1,1,1,CR2,1];

    x1=lsqnonlin(@(s) ss_PCET(s,k1),[x1(1);x1(2);x1(3)],[0 0 0],[1 1 1]);
    if(CR2==0)
        x2=[x1(1) x1(3)];
    end
    
    % Storing reaction rates and coverages in lists
    r1=[r1 k1(7)*x1(1)*x1(3)];
    r2=[r2 k1(8)*x1(2)*x1(3)];
    r3=[r3 k1(7)*x2(1)*x2(2)];
    thetaR1=[thetaR1 x1(1)];
    thetaR2=[thetaR2 x1(2)];
    thetaH=[thetaH x1(3)];    

end

enhance=r1./r3;
select=r1./r2;

figure;
ax1=subplot(2,1,1);
plot(0.1:0.1:4,enhance(2:end));
ylabel('EF')
ax2=subplot(2,1,2);
plot(0:0.1:4,thetaR1,0:0.1:4,thetaR2,0:0.1:4,thetaH);
ylabel('Coverage')

ax1.Position(2)=ax2.Position(2)+ax2.Position(4);
linkaxes([ax1 ax2],'x')


%% Functions

% Parameters
% s - Coverages
% k- Rate parameters

function f=ss_LH3(s,k)
f(1)=k(1)*k(9)*(1-s(1)-s(2)-s(3))-k(2)*s(1)-k(7)*s(1)*s(3);
f(2)=k(3)*k(10)*(1-s(1)-s(2)-s(3))-k(4)*s(2)-k(8)*s(2)*s(3);
f(3)=k(5)*k(11)*(1-s(1)-s(2)-s(3))-k(6)*s(3)-k(8)*s(2)*s(3)-k(7)*s(1)*s(3);
end

function f=ss_ER(s,k)
f(1)=k(1)*k(9)*(1-s(1)-s(2))-k(2)*s(1)-k(7)*s(1)*s(2);
f(2)=k(5)*k(11)*(1-s(1)-s(2))-k(6)*s(2)-k(8)*k(10)*s(2)-k(7)*s(1)*s(2);
end

function f=ss_PCET(s,k) 
f(1)=k(1)*k(9)*(1-s(1)-s(2)-s(3))-k(2)*s(1)-k(7)*s(1)*s(3);
f(2)=k(3)*k(10)*(1-s(1)-s(2)-s(3))-k(4)*s(2)-k(8)*s(2)*k(11);
f(3)=k(5)*k(11)*(1-s(1)-s(2)-s(3))-k(6)*s(3)-k(7)*s(1)*s(3);
end

