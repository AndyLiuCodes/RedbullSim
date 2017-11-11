
function main

clf
g = 6.67*10^(-11);
mEarth = 5.97*10^24;
rEarth = 6378100;  
INI_HEIGHT = 38969.4;
FELIX_MASS = 118;

[time,~,~] = xlsread('data_clean_more_fixed_simplest','data_clean_more_fixed_label','A1:A15348');
[alt,~,~] = xlsread('data_clean_more_fixed_simplest','data_clean_more_fixed_label','B1:B15348');
[v,~,~] = xlsread('data_clean_more_fixed_simplest','data_clean_more_fixed_label','C1:C15348');

speed = v;
%Before parachute data
veloT = max(speed);
increVelo = find(speed == veloT);
heightOfVelo = alt(increVelo(1));
r = 6.38*10^6+heightOfVelo;
newr = r/(6.38*10^6);
gH = 17/(newr^2);
fg = gH*(-9.8*FELIX_MASS);
speed = v/3.6;
maxtime = max(time);


%% Part 1
% Answer some questions here in these comments...
% How accurate is the model for the first portion of the minute? 

%The first portion of the first minute is very accurate until about t = 30
%seconds. At this point, our modeled data vears off from the measured data.
%One reason is since our modeled data does not include air resistance,
%Felix does not hit terminal velocity, so our data will be inaccurate in that section.

% How accurate is the model for the last portion of that first minute? 

%The last portion of the first minute is not as accurate as the first
%portion. As I said above, a reason could be because our modeled data does
%not hit terminal velocity like the measured data because our modeled data
%does not include air resistance.

% Comment on the acceleration calculated from the measured data. 
% Is there any way to smooth the acceleration calculated from the data?

%Our acceleration calculated from the measured data at any point can jump
%between a range of 10ms^-2, meaning that the maximum spike in our graph is
%a 10ms^-2 jump. To lessen this error. one way of smoothing the graph is 
%to gather more data points. With more data points in a given time interval,
%the acceleration data points will be closer to its successor resulting in
%a smoother graph. There is also a smooth function in Matlabs which we took 
%advantage of to clean up our acceleration graph.

part = 1;

plotComparisons(60,'Part 1 - Freefall')

%% Part 2
% Answer some questions here in these comments...
% Estimate your uncertainty in the mass that you have chosen (at the 
%     beginning of the jump). 

% Based on our online research, we chose Felix's total mass to be 118Kg to
% be a good estimate. However, based on our online research, we found a
% good estimate for the minimum range of his mass to be about 100Kg and the
% maximum to be 130Kg. This gives our estimate of 118Kg an error bound of
% about +/- 15Kg

% How sensitive is the velocity and altitude reached after 60 seconds to 
%    changes in the chosen mass?

%After 60 seconds, Felix would have already hit terminal velocity. Terminal
%velocity is dependant on many factors, one of them being mass. Since
%terminal velocity depends on mass, altitude does aswell. However, when
%running our graph with m = 100 instead of m = 118. The changes in the
%velocity and altitude are insignificant and almost unoticable, therefore 
%we see that the graph is not very sensitive to an 18Kg mass change.

part = 2;

plotComparisons(60, 'Part 2 - Simple Air Resistance')
%% Part 3
% Answer some questions here in these comments...
% Felix was wearing a pressure suit and carrying oxygen. Why? 
% What can we say about the density of air in the stratosphere?
% How is the density of air different at around 39,000 meters than it 
% is on the ground?

%Felix wore a full body pressure suit that acted as his life support
%system, providing protection from temperatures of 100 degrees to -90
%degrees F. At very high altitudes, the suit maintains a constant and safe
%pressure from the extreme low pressure of the atmosphere at that height.He
%was carrying oxygen because as altitude increases, the concentration of
%oxygen decreases. The air density at sea level is 1.225kg/m^3 whil the air
%density at the stratosphere and upwards is near 0, therefore the drag due
%to air resistance increases closer to the ground. 

% What are the factors involved in calculating the density of air? 
% How do those factors change when we end up at the ground but start 
% at the stratosphere?  Please explain how calculating air density up 
% to the stratosphere is more complicated than say just in the troposphere.

%The main factors that determine air density is temperature, pressure,
%altitude, and humidity. At sea level, the air density is much greater,
%temperature is warmer and pressure is much greater compared to the stratosphere. The reason calculating
%The air density at the stratosphere is more complicated than calculating the
%density at the troposphere because in the stratosphere, there are 3
%layers to consider the upper stratosphere, lowerstratosphere, and
%stratospause. Which is more complicated than calculating the 
%2 layers in the troposphere (troposphere and tropopause) 


% What method(s) can we employ to estimate [the ACd] product? 

%To estimate the ACd product, we found the ACd product by finding his maximum speed,
% then finding his height and air density at that maximum speed, after we
% derived the formula ACd = 2*Fg/(rho*velocity^2) and solved for ACd
% What is your estimated [ACd] product?
% the result was 0.588

%[Given what we are told in the textbook about the simple drag constant, b,] 
% does the estimate for ACd seem reasonable?

% No, 0.2 is almost 3x less than what the coef. of drag should be for a human in
% freefall. Also, our estimate of 0.588 a an underestimate because we did not
% account for the area of his suit.
part = 3;

plotComparisons(270, 'Part 3 - Fall with Complex Air Resistance')
%% Part 4
% Answer some questions here in these comments...
% What is the actual gravitational field strength around 39,000 meters? 
%   (See Tipler Volume 1 6e page 369.) 
% -9.69ms^-2
% How sensitive is the altitude reached after 4.5 minutes to simpler and 
% more complicated ways of modelling the gravitational field strength? 

%The actual gravitatational field strength at 39,000 meters is
%approximately -9.69ms^-2. Felix was at a height of 2400 meters above the
%Earth's surface 4.5 minutes into the jump. The gravitational field
%strength at this height compared to the strength at the Earth's surface
%differs by less than 0.01ms^-2, therefore implying that the method we use
%to compute the gravitational field strength does not matter and is 
%insignificant because the error is very small.  

% What other changes could we make to our model? Refer to, or at least 
% attempt to explain, the physics behind any changes that you propose. 

%One change we could factor into our model ito make it more accurate is 
%to find out Felix's horizontal component of velocity during the fall. 
%With unpredictable wind patterns etc., Felix would
%not have fallen straight down from his jump position in the atmosphere,
%but instead a considerable distance away from the position directly
%below his platform. Factoring this in would impact our graphs greatly
%because with his horizontal component of velocity the total time to fall
%would have been longer than if he just fell straight down. 


% What is a change that we could make to our model that would result in 
% insignificant changes to the altitude reached after 4.5 minutes? 


%After 4.5 minutes, a change we could make in our model that would make a 
%insignificant change is wether we used just -9.80ms^-2 or used the
%general formula -GMm/r^2 to calculate field strength because 2000 meters
%above the Earth's surface is negligible to the size of Earth's radius.

% How can we decide what change is significant and what change is 
%   insignificant?


%The total distance Felix fell is about 39,000 meters. This implies that
%any change that would be significant to the model would have to change
%significantly from 6.38*10^6 meters(Earth's radius) to 6.38*10^6 + 39,000
%meters(his jump height). Examples are air pressure, air density and
%temperature.



% [What changes did you try out to improve the model?  (Show us your changes
%   even if they didn't make the improvement you hoped for.)]
% <put your answer here in these comments>

%We tried different masses, change in gravity(using the gravitational force
%instead of a constant for g), and different cross sectional areas for Felix.

part = 4;

plotComparisons(270, 'Part 4 - Fall with Air Resistance and Changing Gravity');

%% Part 5
% Answer some questions here in these comments...
% At what altitude does Felix pull the ripcord to deploy his parachute? 

%approximately 2550 meters above Earth's sea level

% Recalculate the ACd product with the parachute open, and modify your 
%   code so that you use one ACd product before and one after this altitude. 
%   According to this version of the model, what is the maximum magnitude 
%   of acceleration that Felix experiences? 

%From our acceleration graph, we observe that the maximum change in
%acceleration that Felix experiences, is about 200ms^-2. 

%   How safe or unsafe would such an acceleration be for Felix?

%Felix experienced a change in acceleration of about 200ms^-2 which is
%about 20G's. To put this in perspective, a lethal car accident ranges from
%30 to 60 G's of Force. From this information, we can see that 20G's is not
%lethal, but also not good for the body. the acceleration was safe enough
%for him to survive.

part = 5;

%Make a single acceleration-plot figure that includes, for each of the 
%model and the acceleration calculated from measurements, the moment when 
%the parachute opens and the following 10 or so seconds. If you have 
%trouble solving this version of the model, just plot the acceleration 
%calculated from measurements. 

SingleAccelPlot

%% Part 6 
% Answer some questions here in these comments...
% How long does it take for Felix’s parachute to open?

%It took about 4.4 seconds for his parachute to open.

part = 6;

%Redraw the acceleration figure from the previous Part but using the new 
%   model. Also, using your plotting function from Part 1, plot the 
%   measured/calculated data and the model for the entire jump from 
%   stratosphere to ground.


plotComparisons(541, 'Part 6 - Entire Redbull Stratos Jump');
SingleAccelPlotRevised
%% nested functions  
% nested functions below are required for the assignment.  
% see Downey Section 10.1 for discussion of nested functions

function res = fall(t, X)
    %FALL <Summary of this function goes here>
    %   <Detailed explanation goes here>

    % do not modify this function unless required by you for some reason! 

    p = X(1); % the first element is position
    v = X(2); % the second element is velocity

    dpdt = v; % velocity: the derivative of position w.r.t. time
    dvdt = acceleration(t, p, v); % acceleration: the derivative of velocity w.r.t. time

    res = [dpdt; dvdt]; % pack the results in a column vector
end

function res = acceleration(t, p, v)
    % <insert description of function here>
    % input...
    % t: time
    % p: position
    % v: velocity
    % output...
    % res: acceleration

    % do not modify this function unless required by you for some reason! 

    a_grav = gravityEst(p); 

    if part == 1 % variable part is from workspace of function main.
        res = -a_grav;
    else
        m = mass(t, v);
        b = drag(t, p, v, m);

        f_drag = b * v^2;
        a_drag = f_drag / m;
        res = -a_grav + a_drag;
    end
end

% Please paste in or type in code into the below functions as may be needed.

function a_grav = gravityEst(p)
    % estimate the acceleration due to gravity as a function of altitude, p
    A_GRAV_SEA = 9.807;  % acceleration of gravity at sea level in m/s^2

     if part <= 3
        a_grav = A_GRAV_SEA;
     else
        a_grav = g*(mEarth)/(p + rEarth)^2;
     end
     
end

function res = mass(t, v)
    
    m = FELIX_MASS;
    
    res = m;
end

function res = drag(t, p, v, m)
% <insert description of function here>
     
    if part == 2
        res = 0.2;
    else
        rho = stdatmo(heightOfVelo);
        ACd = 2*abs(fg)/(rho*veloT^2);
        rhoC = stdatmo(p);

        
        if t >= 264 && part > 4
        Cd = ACd /(0.7*1.7-.1);
        ACd = 25*Cd;  %Area of parachute * Cd  
        end
        
        res = .5*rhoC*ACd;
        
        
    end
end

%% Additional nested functions
% Nest any other functions below.  
%Do not put functions in other files when you submit.

function output = SingleAccelPlot

lengthTime = length(time);

for j=2:(lengthTime)
   actAccel(j) = (-1)*(speed(j)-speed(j-1))/(time(j)-time(j-1));
end

        
        
[t,a] = ode45(@fall,[0,270],[38969.4,0]);
        

fspeed = a(:,2);

%modelled acceleration

lengthTmod = length(t);

for k=2:(lengthTmod )
   modAccel(k) =  (fspeed(k)-fspeed(k-1))/((t(k) - t(k-1)));
end
        figure(4);
        clf
        plot(time,actAccel,'ro-');
        hold on 
        plot(t,modAccel,'bx-');
        title('Part 5 - Single Acceleration Plot')
        ylabel('Acceleration(m/s^2)');
        xlabel('Time(s)');
        legend('Measured Data', 'Modeled Data','Location','Best');
        xlim([250,274])
        output = gcf;
end


function output = SingleAccelPlotRevised

lengthTime = length(time);

for j=2:(lengthTime)
   actAccel(j) = (-1)*(speed(j)-speed(j-1))/(time(j)-time(j-1));
end

        
        
[t,a] = ode45(@fall,[0,270],[38969.4,0]);
        

fspeed = a(:,2);


%modelled acceleration

lengthTmod = length(t);

for k=2:(lengthTmod )
   modAccel(k) =  (fspeed(k)-fspeed(k-1))/((t(k) - t(k-1)));
end
        figure(4);
        clf
        plot(time,actAccel,'ro-');
        hold on 
        plot(t,modAccel,'bx-');
        title('Part 6 - Single Acceleration Plot Revised')
        ylabel('Acceleration(m/s^2)');
        xlabel('Time(s)');
        legend('Measured Data', 'Modeled Data','Location','Best');
        xlim([260,264.4])
        ylim([-36,60]);
        output = gcf;
end


function output = plotComparisons(max, PartTitle)
        
hold off  
%actual acceleration
lengthTime = length(time);
for j=2:(lengthTime)
   actAccel(j) = (-1)*(speed(j)-speed(j-1))/(time(j)-time(j-1));
end
%modeled

[t,a] = ode45(@fall,[0,max],[INI_HEIGHT,0]);
falt = a(:,1);
fspeed = a(:,2);
hold on 
%modelled acceleration
lengthTmod = length(t);
for k=2:(lengthTmod )
   modAccel(k) =  (fspeed(k)-fspeed(k-1))/((t(k) - t(k-1)));
end

%Altitude
figure(1)
clf
hold off
plot(time,alt)
hold on 
plot(t,falt)
title(PartTitle)
ylabel('Altitude(m)')
xlabel('Time(s)');
legend('Measured Data','Modeled Data','Location','Best');
if part <= 2
xlim([0,60])
elseif part > 2 && part <=4
xlim([0,270])
else
xlim([0,maxtime])
end


%Velocity
figure(2)
clf
hold off
plot(time,speed)
hold on 
plot(t,abs(fspeed));  
title(PartTitle)
xlabel('Time(s)');
ylabel('Velocity(m/s)');
legend('Measured Data', 'Modeled Data','Location','Best');
if part <= 2
xlim([0,60])
elseif part > 2 && part <=4
xlim([0,270])
else
xlim([0,maxtime])
end


%Acceleration
figure(3)
clf
hold off
plot(time,smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(actAccel)))))))))));
hold on 
if part == 1
plot(t,ones(length(t),1)*(-9.8))
else
plot(t,modAccel)
end
title(PartTitle)
ylabel('Acceleration(m/s^2)');
xlabel('Time(s)');
legend('Measured Data', 'Modeled Data','Location','Best');
if part <= 2
xlim([0,60])
elseif part > 2 && part <=4
xlim([0,270])
else
xlim([0,maxtime])
end
if part == 6
ylim([-25,25])
end
        
        
output = gcf;
        
end


% end of nested functions
end % closes function main.  
