
%% Matlab set up

clear all
clc

Width_array  = input('Amount of horizontal (leftside) lasers: ');
Height_array = input('Amount of vertical (upperhalf) lasers: ');

%% Matrix Sizes

book=26.6; %(cell size=13.3)
matrix_length=100; %
matrix_width=62; % SOS!! to be added
wall_height=37; %SOS change

%% Obeject positions
box_y=29.5; %wooden box
box_x=26;

cards_y=8.4; %Poker cards
cards_x=48.9;

card_boxes_y=8.8; %Perudo game
card_boxes_x=73.4;

laptop_y=32.3; %Laptop keyboard
laptop_x=70.1;

laptop_screen_y=56.5;
laptop_screen_x=70.1;

camera_y=50.5;
camera_x=71.5;
camera_h=24;


%% Object sizes

box_s_x=10.4;
box_s_y=10.4;
box_s_h=15.5;

laptop_s_x=35.5;
laptop_s_y=24.3;
laptop_s_h=1.3;

laptop_screen_s_x=35.5;
laptop_screen_s_y=0.4;
laptop_screen_s_h=24.5;

cards_s_x=6.1;
cards_s_y=8.9;
cards_s_h=7.8;

wood_s_x=25;
wood_s_y=18.3;
wood_s_h=13;

card_boxes_s_x=10.4;
card_boxes_s_y=10.4;
card_boxes_s_z=15.3;

%% Matrix object implementations

scale=10;% we add a scale to get ride of milimeters to create proper matlab cells

global_matrix_x = matrix_length*scale;
global_matrix_y = matrix_width*scale;

global_matrix (1:global_matrix_y,1:global_matrix_x)= 1;


for c= box_x*scale:1:(box_x+box_s_x)*scale
    for i= global_matrix_y- ((box_y+box_s_y)*scale):1:global_matrix_y-(box_y*scale)
        global_matrix(i,c)=scale*box_s_h;  %We assign ones to the cells where the box belongs in the matrix
    end
end


for c= cards_x*scale:1:(cards_x+cards_s_x)*scale
    for i= global_matrix_y- ((cards_y+cards_s_y)*scale):1:global_matrix_y-(cards_y*scale)
        global_matrix(i,c)=scale*cards_s_h;  %We assign ones to the cells where the box belongs in the matrix
    end
end


for c= round(card_boxes_x*scale):1:round((card_boxes_x+card_boxes_s_x)*scale)
    for i= global_matrix_y- (round((card_boxes_y+card_boxes_s_y)*scale)):1:global_matrix_y-round(card_boxes_y*scale)
        global_matrix(i,c)=round(scale*card_boxes_s_z);  %We assign ones to the cells where the box belongs in the matrix
    end
end

for c= round(laptop_x*scale):1:round((laptop_x+laptop_s_x)*scale)
    for i= global_matrix_y- (round((laptop_y+laptop_s_y)*scale)):1:global_matrix_y-round(laptop_y*scale)
        global_matrix(i,c)=round(scale*laptop_s_h);  %We assign ones to the cells where the box belongs in the matrix
    end
end


for c= round(laptop_screen_x*scale):1:round((laptop_screen_x+laptop_screen_s_x)*scale)
    for i= global_matrix_y- (round((laptop_screen_y+laptop_screen_s_y)*scale)):1:global_matrix_y-round(laptop_screen_y*scale)
        global_matrix(i,c)=round(scale*laptop_screen_s_h);  %We assign ones to the cells where the box belongs in the matrix
    end
end

%% Laser Beam

camera_height = scale* camera_h;
Laser_horizontal_x=scale*matrix_length;
Laser_horizontal_y= scale*(camera_y+matrix_width);
Laser_field_angle_x= atan2(Laser_horizontal_x/2,Laser_horizontal_y);

Laser_vertical_h=scale*48;%Height from camera, to wall highest part recorded by camera (remember that the vertical height from table to max recorded height is 72 cm)
Laser_vertical_y= scale*(camera_y+matrix_width);
Laser_field_angle_y= atan2(Laser_vertical_h,Laser_vertical_y);

global_laser_matrix(1:scale*matrix_width,1:scale*matrix_length)= NaN;
upper_matrices(1:scale*matrix_width,1:scale*matrix_length)= NaN;
lower_matrices(1:scale*matrix_width,1:scale*matrix_length)= NaN;


[mm ,nn] = size(global_matrix);
numb_laser_beams_per_side=Width_array; %How many lasers are on the right side
numb_laser_beams_per_upper_height=Height_array; %How many lasers are on the upper side

Local_laser_x_angle= Laser_field_angle_x/numb_laser_beams_per_side ;
Local_laser_y_angle= Laser_field_angle_y/numb_laser_beams_per_upper_height;

% first_detection_boulean(1:mm,1:nn)=true;
first_detection_boulean=true(mm,nn,4); % 1 is for left-up, 2 is for right-up, 3 is for left-down, 4 is for right-down

contour_matrix(1:numb_laser_beams_per_upper_height*2,1:numb_laser_beams_per_side*2)=NaN;


for a=1:numb_laser_beams_per_upper_height


    for i= 1:1:mm



        for n=1:numb_laser_beams_per_side %Repeats with each laser beam

            Y(i,n)= (i+ (scale*(camera_y)) ) /(cos(Local_laser_x_angle*n)); % Y is distance that the light has travelled (Not perpendicular!!)
            X(i,n)= Y(i,n)*sin(Local_laser_x_angle*n); % This is X direction (lengthwise) that the light hits for each row of the global matrix



            Left_laser(i,n)=round(1+(Laser_horizontal_x/2)-X(i,n)); %X ccordinate of cell in table
            Right_laser(i,n)=round((Laser_horizontal_x/2)+X(i,n));  %Y coordinate of cell in table

            Z_Y(i,n)= (i+ (scale*(camera_y)) ) /(cos(Local_laser_y_angle*a)); %Hypothenus lenght (upwardslaserlinelength)
            Z_X(i,n)=  Z_Y(i,n)*sin(Local_laser_y_angle*a);




            upper_matrices(1+mm-i,Left_laser(i,n) ) =round(camera_height + Z_X(i,n)) ; % the value is its height
            upper_matrices(1+mm-i,Right_laser(i,n) ) = round(camera_height + Z_X(i,n)) ; % the value is its height

            lower_matrices(1+mm-i,Left_laser(i,n) )=round(camera_height -Z_X(i,n)); % the value is its height
            lower_matrices(1+mm-i,Right_laser(i,n) ) = round(camera_height -Z_X(i,n)); % the value is its height



            global_laser_matrix(1+mm-i,Left_laser(i,n))=2;
            global_laser_matrix(1+mm-i,Right_laser(i,n))=2;




            if first_detection_boulean(a,n,1)==true

                if ~isnan(global_matrix(1+mm-i,Left_laser(i,n)) ) && upper_matrices(1+mm-i,Left_laser(i,n)) >= 1 &&  upper_matrices(1+mm-i,Left_laser(i,n)) <= global_matrix(1+mm-i,Left_laser(i,n))%upper_matrices is the height of the laser %max Height object value of where lidar strikes in global matrix (x,y)
                    % THe first statement finds if global matrix position where lidar strikes
                    % are not NaN (x,y), the second statement checks if the laser is higher
                    % than table height, the thrid statement checks if laser strikes the top
                    % upper part of the object and lower


                    %contour: x,y = depth

                    contour_matrix(numb_laser_beams_per_upper_height+a,1+numb_laser_beams_per_side-n)=1+mm-i;

                    first_detection_boulean(a,n,1)=false;

                end

            end


            if first_detection_boulean(a,n,2)==true
                if ~isnan(global_matrix(1+mm-i,Right_laser(i,n)) ) && upper_matrices(1+mm-i,Right_laser(i,n) )  >= 1 &&  upper_matrices(1+mm-i,Right_laser(i,n)) <= global_matrix(1+mm-i,Right_laser(i,n))

                    contour_matrix(numb_laser_beams_per_upper_height+a,numb_laser_beams_per_side+n)=1+mm-i;

                    first_detection_boulean(a,n,2)=false;

                end
            end

            if first_detection_boulean(a,n,3)==true
                if ~isnan(global_matrix(1+mm-i,Left_laser(i,n)) ) && lower_matrices(1+mm-i,Left_laser(i,n))   >= 1 &&  lower_matrices(1+mm-i,Left_laser(i,n)) <= global_matrix(1+mm-i,Left_laser(i,n))


                    contour_matrix(1+numb_laser_beams_per_upper_height-a,1+numb_laser_beams_per_side-n)=1+mm-i;
                    first_detection_boulean(a,n,3)=false;
                end
            end


            if first_detection_boulean(a,n,4)==true
                if ~isnan(global_matrix(1+mm-i,Right_laser(i,n)) )  && lower_matrices(1+mm-i,Right_laser(i,n) )  >= 1 &&  lower_matrices(1+mm-i,Right_laser(i,n)) <= global_matrix(1+mm-i,Right_laser(i,n))

                    contour_matrix(1+numb_laser_beams_per_upper_height-a,numb_laser_beams_per_side+n)=1+mm-i;
                    first_detection_boulean(a,n,4)=false;
                end
            end




        end



    end




    depth_up_matrices{a}=upper_matrices;
    depth_low_matrices{a}=lower_matrices;
end

%% Contour plot

figure(2); clf(2);
[xi, yi] = meshgrid(1:size(contour_matrix, 2), 1:size(contour_matrix, 1));
colormap(hot)
contourf(xi,yi,contour_matrix,100)
%colormap(gray);
colorbar;
axis equal
