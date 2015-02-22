function draw_ellipse(x, y, a, b, angle, a1, a2,w)%, steps) 


%if (steps == 0)
    steps = 36;
%end

if a1<0, a1 = a1+2*pi; end
if a2<0, a2 = a2+2*pi; end
if a2<a1, a2 = a2+2*pi; end
a1 = a1*180/pi;
a2 = a2*180/pi;
points = [];
% angle = -20;
% x = 5;
% y = 4;
% a = 3;
% b = 1.4;
  
%beta = angle * (pi / 180); 
beta = angle;
sinbeta = sin(beta);
cosbeta = cos(beta);


%for i = 0:360/steps:360 
for i = a1:(a2-a1)/steps:a2 
    alpha = i * (pi / 180) ;
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);
 
    X = x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta);
    Y = y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta);
    points = [points; X Y];
end
% alpha1 = [0; pi];
% sinalpha1 = sin(alpha1);
% cosalpha1 = cos(alpha1);
% X = x + (a * cosalpha1 * cosbeta - b * sinalpha1 * sinbeta);
% Y = y + (a * cosalpha1 * sinbeta + b * sinalpha1 * cosbeta);
% plot(X, Y,'*r');
hold on
if a==b 
   plot(points(:,1), points(:,2),'r','LineWidth',w);
else
    plot(points(:,1), points(:,2),'m','LineWidth',w);
end
    
end