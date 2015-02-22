function display_polyrect( poly, pts )
%DISPLAY_POLY Summary of this function goes here
%   Detailed explanation goes here

%imshow(im), hold on;

if ~isempty(pts)
    plot(pts(:,1)+1, pts(:,2)+1,'*r');
end

for i=1:size(poly,1)
    display_rect(poly(i,:));    
end

end


function display_rect(r)

x1 = r(1) - r(6) * r(7);
y1 = r(2) + r(5) * r(7);
x2 = r(3) - r(6) * r(7);
y2 = r(4) + r(5) * r(7);
x3 = r(3) - r(6) * r(8);
y3 = r(4) + r(5) * r(8);
x4 = r(1) - r(6) * r(8);
y4 = r(2) + r(5) * r(8);

text(x1+1,y1+1,'1');
text(x2+1,y2+1,'2');
text(x3+1,y3+1,'3');
text(x4+1,y4+1,'4');
line([x1+1 x2+1],[y1+1 y2+1],'Color','g');
line([x2+1 x3+1],[y2+1 y3+1],'Color','g');
line([x3+1 x4+1],[y3+1 y4+1],'Color','g');
line([x4+1 x1+1],[y4+1 y1+1],'Color','g');
end

