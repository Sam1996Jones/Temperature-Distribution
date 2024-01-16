clear;
clc;
Delta_x = 5; % change in x direction [cm]
Delta_y = 5; % change in y direction [cm]
Delta_t = 10; % time step size [s]
Diff = 0.2; % Thermal diffusivity [cm^2/s]
%T = (0:Delta_t:600);
Y = zeros(6,8); % 48 nodes 
Plate = zeros(6,8); % Size of plate in a point array (y,x)
Plate(1:6,1:8) = 40; % inital temperature of plate 40 degrees
Plate(:,1) = 50; % set boundary temperature 50 degrees
Plate(6,2:8) = 45; % set boundary temperature 45 degrees 
Plate(1:6,8) = 30; % set boundary temperautre 30 degrees 
%Plate(4:5,4:5) = 40; % set boundary temperature 40 degrees
Plate0 = Plate; %Plate0 is previous time step e.g n=0. Plate is current time step e.g n = 0.5 etc

pos = 2; % postion loop 
pos1(1)=Plate(5,2); %Position vector for certain points on plate (in question)
pos2(1)=Plate(2,3);
pos3(1)=Plate(4,6);
pos4(1)=Plate(3,5);
pos5(1)=Plate(2,5);
pos6(1)=Plate(1,5);

Lambda = 2*((Diff*Delta_t)/(Delta_x^2));
a = -Lambda/4; % ADI constants-delta x and delta y are equal
b = 1 + (Lambda/2);
c = -Lambda/4;
A = a;
B = b;
C = c;

for time = Delta_t:Delta_t:600  %Length of simulation in seconds

%%
%%%%%%%%%%%%%%%%% n : 0 -> 0.5 etc
%% Calculate Temperates at boundary D using central difference (Nueman Boundary)
% j = 1
syms v2 v3 v4 v5 v6 v7
for i = 2:7
    p(i) = (Plate0(1,i)*(1-Lambda/2)) + (Lambda/2)*Plate0(2,i) - (Lambda/4); %constant
end
eqn2 = v2*b + v3*c == p(2) - Plate(1,1)*a; %simultaneous eqns
eqn3 = v2*a + v3*b + v4*c == p(3);
eqn4 = v3*a + v4*b + v5*c == p(4);
eqn5 = v4*a + v5*b + v6*c == p(5);
eqn6 = v5*a + v6*b + v7*c == p(6);
eqn7 = v6*a + v7*b == p(3) - Plate(1,8)*c;
sol = solve([eqn2 eqn3 eqn4 eqn5 eqn6 eqn7], [v2 v3 v4 v5 v6 v7]);
Plate(1,2) = sol.v2;
Plate(1,3) = sol.v3;
Plate(1,4) = sol.v4;
Plate(1,5) = sol.v5;
Plate(1,6) = sol.v6;
Plate(1,7) = sol.v7;

%% j = 2
for i = 2:7
    p(i) = (Plate0(2,i)*(1-Lambda/2)) + (Lambda/4)*(Plate0(1,i)+Plate0(3,i)); % constant
end
eqn2 = v2*b + v3*c == p(2) - Plate(2,1)*a; %simultaneous eqns
eqn3 = v2*a + v3*b + v4*c == p(3);
eqn4 = v3*a + v4*b + v5*c == p(4);
eqn5 = v4*a + v5*b + v6*c == p(5);
eqn6 = v5*a + v6*b + v7*c == p(6);
eqn7 = v6*a + v7*b == p(3) - Plate(2,8)*c;
sol = solve([eqn2 eqn3 eqn4 eqn5 eqn6 eqn7], [v2 v3 v4 v5 v6 v7]);
Plate(2,2) = sol.v2;
Plate(2,3) = sol.v3;
Plate(2,4) = sol.v4;
Plate(2,5) = sol.v5;
Plate(2,6) = sol.v6;
Plate(2,7) = sol.v7;
%% j = 3, including boundary F
for i = 2:3
    p(i) = (Plate0(3,i)*(1-Lambda/2)) + (Lambda/4)*(Plate0(2,i)+Plate0(4,i)); % constant
end
for i = 4:5
    p(i) = (Plate0(4,i)*(1-Lambda/2)) + (Lambda/2)*(Plate0(2,i)-1); % constant
end
for i = 6:7
    p(i) = (Plate0(3,i)*(1-Lambda/2)) + (Lambda/4)*(Plate0(2,i)+Plate0(4,i)); % constant
end
eqn2 = v2*b + v3*c == p(2) - Plate(3,1)*a; %simultaneous eqns
eqn3 = v2*a + v3*b + v4*c == p(3);
eqn4 = v3*a + v4*b + v5*c == p(4);
eqn5 = v4*a + v5*b + v6*c == p(5);
eqn6 = v5*a + v6*b + v7*c == p(6);
eqn7 = v6*a + v7*b == p(3) - Plate(3,8)*c;
sol = solve([eqn2 eqn3 eqn4 eqn5 eqn6 eqn7], [v2 v3 v4 v5 v6 v7]);
Plate(3,2) = sol.v2;
Plate(3,3) = sol.v3;
Plate(3,4) = sol.v4;
Plate(3,5) = sol.v5;
Plate(3,6) = sol.v6;
Plate(3,7) = sol.v7;

%%
syms x y
for j = 4:5 %rows 4 -> 5
   for i = 2:3 %left side of hole
       d(i) = Plate0(j,i)*(1-Lambda/2) + (Lambda/4)*(Plate0(j-1,i)+Plate0(j+1,i));
   end
   eqn1 = b*x + c*y == d(2) - a*Plate(j,1);
   eqn2 = a*x + b*y == d(3) - c*Plate(j,4);
   sol2 = solve([eqn1 eqn2], [x y]);
   Plate(j,2) = sol2.x;
   Plate(j,3) = sol2.y;
   for i = 6:7 %right side of hole
       d(i) = Plate0(j,i)*(1-Lambda/2) + (Lambda/4)*(Plate0(j-1,i)+Plate0(j+1,i));
   end
   eqn1 = b*x + c*y == d(6) - a*Plate(j,5);
   eqn2 = a*x + b*y == d(7) - c*Plate(j,8);
   sol2 = solve([eqn1 eqn2], [x y]);
   Plate(j,6) = sol2.x;
   Plate(j,7) = sol2.y;
end
Plate0 = Plate; % Plate0 n:0.5 and Plate becomes n:1 etc

%%
%%%%%%%%%%%%%%%%%%% n : 0.5 -> 1 etc

% i = 2

for i = 2:3
    for j = 2:5
       D(j) = Plate0(j,i)*(1-Lambda/2) + (Lambda/4)*(Plate0(j,i-1)+Plate0(j,i+1)); 
    end
    eqn2 = B*v2 + c*v3 == D(2) - A*Plate0(1,i);
    eqn3 = A*v2 + B*v3 + C*v4 == D(3);
    eqn4 = A*v3 + B*v4 + C*v5 == D(4);
    eqn5 = A*v4 + B*v5 == D(5) - C*Plate(6,i);
    sol = solve([eqn2 eqn3 eqn4 eqn5], [v2 v3 v4 v5]);
    Plate(2,i) = sol.v2;
    Plate(3,i) = sol.v3;
    Plate(4,i) = sol.v4;
    Plate(5,i) = sol.v5;
    %boundary D
    D(1) = Plate0(1,i)*(1-Lambda/2) + (Lambda/4)*(Plate0(1,i-1)+Plate0(1,i+1));
    Plate(1,i) = ( D(1)- C*Plate(2,2) - A*(Plate(2,2)-1) )/B; 
    %
end
for i = 6:7
    for j = 2:5
       D(j) = Plate0(j,i)*(1-Lambda/2) + (Lambda/4)*(Plate0(j,i-1)+Plate0(j,i+1)); 
    end
    eqn2 = B*v2 + c*v3 == D(2) - A*Plate0(1,i);
    eqn3 = A*v2 + B*v3 + C*v4 == D(3);
    eqn4 = A*v3 + B*v4 + C*v5 == D(4);
    eqn5 = A*v4 + B*v5 == D(5) - C*Plate(6,i);
    sol = solve([eqn2 eqn3 eqn4 eqn5], [v2 v3 v4 v5]);
    Plate(2,i) = sol.v2;
    Plate(3,i) = sol.v3;
    Plate(4,i) = sol.v4;
    Plate(5,i) = sol.v5;
    %Boundary D
    D(1) = Plate0(1,i)*(1-Lambda/2) + (Lambda/4)*(Plate0(1,i-1)+Plate0(1,i+1));
    Plate(1,i) = ( D(1)- C*Plate(2,2) - A*(Plate(2,2)-1) )/B;
    %
end

syms z
% between boundary D and F
for i = 4:5
    for j = 1:3
        D(j) = Plate0(j,i) + (Lambda/4)*(Plate0(j,i-1) - 2*Plate0(j,i) + Plate0(j,i+1));
    end
    eqn1 = (A+C)*y + B*x == D(1)+A;
    eqn2 = A*x + B*y + C*z == D(2);
    eqn3 = (A+C)*y + B*z == D(3) + 2*C;
    sol = solve([eqn1 eqn2 eqn3], [x y z]);
    Plate(1,i) = sol.x;
    Plate(2,i) = sol.y;
    Plate(3,i) = sol.z;
end

Plate0 = Plate;

pos1(pos)=Plate(5,2); %Position vector for certain points on plate (in question)
pos2(pos)=Plate(2,3);
pos3(pos)=Plate(4,6);
pos4(pos)=Plate(3,5);
pos5(pos)=Plate(2,5);
pos6(pos)=Plate(1,5);
pos = pos + 1;
 
disp(time)  %displays time step

if time == 20 || time == 90 || time == 300 
    figure();
    contourf(0:5:35,0:5:25,flip(Plate)); %contour plot, flip maxtrix to allow contour to look like figure 1
    colbar = colorbar;
    ylabel(colbar, 'Temperature [\circ]')
    str = sprintf('Contour Plot at %ds', time);
    title(str);
    xlabel('Width along plate [cm]');
    ylabel('Length along plate [cm]');
    hold on;
    rectangle('Position',[15,5,5,10],'FaceColor','white');
    hold on;
end
end


figure();
contourf(0:5:35,0:5:25,flip(Plate)); %contour plot, flip maxtrix to allow contour to look like figure 1
colbar = colorbar;
ylabel(colbar, 'Temperature [\circ]')
str = sprintf('Contour Plot at %ds', time);
title(str);
xlabel('Width along plate [cm]');
ylabel('Length along plate [cm]');
hold on;
rectangle('Position',[15,5,5,10],'FaceColor','white')
hold on;

figure();
x = 0:10:time;
plot(x,pos1,'y');
xlabel('Time [s]');
ylabel('Temperature [\circ]');
title('Temperature vs Time');
Leg = legend({'1'},'Location','southeast');
title(Leg, 'Position');

figure();
x = 0:10:time;
plot(x,pos2,'m');
xlabel('Time [s]');
ylabel('Temperature [\circ]');
title('Temperature vs Time');
Leg = legend({'2'},'Location','southeast');
title(Leg, 'Position');

figure();
x = 0:10:time;
plot(x,pos3,'c');
xlabel('Time [s]');
ylabel('Temperature [\circ]');
title('Temperature vs Time');
Leg = legend({'3'},'Location','southeast');
title(Leg, 'Position');

figure();
x = 0:10:time;
plot(x,pos4,'r');
xlabel('Time [s]');
ylabel('Temperature [\circ]');
title('Temperature vs Time');
Leg = legend({'4'},'Location','southeast');
title(Leg, 'Position');

figure();
x = 0:10:time;
plot(x,pos5,'g');
xlabel('Time [s]');
ylabel('Temperature [\circ]');
title('Temperature vs Time');
Leg = legend({'5'},'Location','southeast');
title(Leg, 'Position');

figure();
x = 0:10:time;
plot(x,pos6,'b');
xlabel('Time [s]');
ylabel('Temperature [\circ]');
title('Temperature vs Time');
Leg = legend({'6'},'Location','southeast');
title(Leg, 'Position');