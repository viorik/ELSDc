function parse_tmp( )
%PARSE_TMP Summary of this function goes here
%   Detailed explanation goes here
close all
infile = 'tmp.txt';
f = fopen(infile,'r');
im = imread('imtest/cercles.pgm');
imshow(im), hold on;
while ~feof(f)
    s1 = fgetl(f);
    if strcmp(s1,'refined_rectangle')
        r = fscanf(f,'%f %f %f %f %f %f %f %f\n',8);
        display_polyrect(r',[]);
        ok=1;
    elseif strcmp(s1,'CC') || strcmp(s1,'EE')
        ce = fscanf(f,'%f %f %f %f %f \n',5);
        if ce(3)<=0 || ce(4)<=0
            disp ('ellipse not valid')
        else
            draw_ellipse(ce(1)+1,ce(2)+1,ce(3),ce(4),ce(5),0,2*pi,1.3);
        end
        ok=1;
    elseif strcmp(s1,'VALID!!!!!')
        typ = fscanf(f,'%c\n',1);
        if typ=='c'
            ce = fscanf(f,'%f %f %f %f %f %f %f\n',7);
            draw_ellipse(ce(1)+1,ce(2)+1,ce(3),ce(4),ce(5),ce(6),ce(7),2);
        else
            dim = fscanf(f,'%d\n',1);
            pp = fscanf(f,'%f ', dim*2);
            display_poly(pp);
        end
            
    end    
end
fclose(f);

end

