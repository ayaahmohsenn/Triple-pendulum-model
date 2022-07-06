t=linspace(-5,5,100);

for i=1:length(t)
y= -0.18664*cos(2*t-pi/2);
plot(t,y);
comet(t,y)
pause(0.1)
   end