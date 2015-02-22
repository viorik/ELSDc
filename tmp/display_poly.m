function display_poly( poly )
%DISPLAY_POLY Summary of this function goes here
%   Detailed explanation goes here

%imshow(im), hold on;

for i=1:4:size(poly,1)
    line([poly(i) poly(i+2)],[poly(i+1) poly(i+3)],'Color','g','LineWidth',2);    
end

end
