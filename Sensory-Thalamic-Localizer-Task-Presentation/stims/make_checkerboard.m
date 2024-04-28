clear pi
hival = 1;
loval = 0;

numchecks = 16;
numbands = 12;
lth = 1081;


intrvl = (360/(numchecks*4));
rintrvl = (lth-1)/numbands;


IMG = zeros(lth) + .5;



lo = true;
count = 1;
for x = 1:length(IMG)
    for  y = 1:length(IMG)
        leng = sqrt((x-1)^2 + (y-1)^2);
        ang = atan((x-1)/(y-1))*(180/pi);
        IND1 = mod(floor(ang/intrvl),2);
        IND2 = mod(floor(leng/rintrvl),2);
        
        if floor(leng/rintrvl)==0 || floor(leng/rintrvl) > numbands - 1
            IMG(y,x) = .5;
        elseif (IND1 && IND2) || (~IND1 && ~IND2)
            IMG(y,x) = loval;
        else
            IMG(y,x) = hival;
        end
    end
end


fullimg = [1 - IMG(end:-1:2,:); IMG];
fullimg = [1 - fullimg(:,end:-1:2) fullimg];

imwrite(fullimg,'checkerboard.tiff');
imwrite(1-fullimg,'checkerboardInv.tiff');