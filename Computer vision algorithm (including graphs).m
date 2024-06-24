clear all;
clc;
close all;


folder_paths = {
    'C:\Users\User\Documents\Thesis\new_colors1', ...
    'C:\Users\User\Documents\Thesis\new_colors2', ...
    'C:\Users\User\Documents\Thesis\new_colors3', ...
    'C:\Users\User\Documents\Thesis\new_colors4', ...
    'C:\Users\User\Documents\Thesis\new_colors5', ...
    'C:\Users\User\Documents\Thesis\new_colors6'
};

% distance matrix 
distance_matrix = zeros(17, 6); 

% Loop over all folder paths
for path_index = 1:numel(folder_paths)
    folder_path = folder_paths{path_index};
    file_list = dir(fullfile(folder_path, '*.png'));
    file_numbers = cellfun(@(x) sscanf(x, 'sim_%d.jpg'), {file_list.name});
    [~, sorted_indices1] = sort(file_numbers);

    n = 0;

    for i = 1:numel(file_list)/2
        tic;
        sorted_indices2 = sorted_indices1(1, n*2 + 1);

        n = i;

        sorted_indices3 = sorted_indices1(1, i*2);

        filename1 = fullfile(folder_path, file_list(sorted_indices2).name);
        filename2 = fullfile(folder_path, file_list(sorted_indices3).name);

        sim1 = imread(filename1); 
        sim2 = imread(filename2); 
        scale_factor = 1;
        sim1_resized = imresize(sim1, scale_factor);
        sim2_resized = imresize(sim2, scale_factor);

        % GrayScale Convertion
        sim1_gray = rgb2gray(sim1_resized);
        sim2_gray = rgb2gray(sim2_resized);

        % Detect corners-> haris
        corners1 = detectHarrisFeatures(sim1_gray, 'MinQuality', 0.1); % Corners in sim1
        corners2 = detectHarrisFeatures(sim2_gray, 'MinQuality', 0.1); % Corners in sim2
        [features1, valid_points1] = extractFeatures(sim1_gray, corners1);
        [features2, valid_points2] = extractFeatures(sim2_gray, corners2);
        indexPairs = matchFeatures(features1, features2);
        matchedPoints1 = valid_points1(indexPairs(:, 1), :);
        matchedPoints2 = valid_points2(indexPairs(:, 2), :);

        % distance calculation
        distances = sqrt(sum((matchedPoints1.Location - matchedPoints2.Location).^2, 2));

        % tresold
        movement_threshold = 10; 
        moved_indices = distances > movement_threshold;
        movedPoints1 = matchedPoints1(moved_indices);
        movedPoints2 = matchedPoints2(moved_indices);
        figure(1); clf(1);
        showMatchedFeatures(sim1_resized, sim2_resized, movedPoints1, movedPoints2, 'montage');
        title('Moved  points');

        mean_x_image_1 = sum(movedPoints1.Location(:,1)) / length(movedPoints1.Location(:,1));
        mean_y_image_1 = sum(movedPoints1.Location(:,2)) / length(movedPoints1.Location(:,2));
        mean_x_image_2 = sum(movedPoints2.Location(:,1)) / length(movedPoints1.Location(:,1));
        mean_y_image_2 = sum(movedPoints2.Location(:,2)) / length(movedPoints1.Location(:,2));

        %mean avg coordinates calculation
        distance_pixels = sqrt((mean_x_image_2 - mean_x_image_1)^2 + (mean_y_image_2 - mean_y_image_1)^2);

        photosize = size(sim1_gray);
        distance_pixels_ratio = distance_pixels / photosize(1,2); %pixels_per_cm;
        lidar_matrix(i, path_index) = 100 * distance_pixels_ratio;

        elapsed_time_lidar(i, path_index) = toc;
    end
end

% Plot the distance matrix
figure(2); clf(2);
plot(1:17, lidar_matrix, '-o');
xlabel('Frames');
ylabel('Distance moved(cm)');
legend('Path A1', 'Path B1', 'Path C1', 'Path D1', 'Path E1', 'Path F1');
title('Distance of object moved ');
grid on;


% Define the single folder path
folder_path = 'C:\Users\User\Documents\Thesis\LargeBook';

scale_factors = [0.3, 0.6, 0.9, 1.2, 1.5, 1.8];

distance_matrix = zeros(17, numel(scale_factors)); 
file_list = dir(fullfile(folder_path, '*.jpeg'));
file_numbers = cellfun(@(x) sscanf(x, 'sim%d.jpg'), {file_list.name});
[~, sorted_indices] = sort(file_numbers);

% Loop over scale factors
for scale_index = 1:numel(scale_factors)
    scale_factor = scale_factors(scale_index);

    for i = 1:numel(file_list)/2
        tic;

        % Get the sorted indices
        idx1 = sorted_indices(2*i-1);
        idx2 = sorted_indices(2*i);

        filename1 = fullfile(folder_path, file_list(idx1).name);
        filename2 = fullfile(folder_path, file_list(idx2).name);

        sim1 = imread(filename1); % Read images
        sim2 = imread(filename2); 
        sim1_resized = imresize(sim1, scale_factor);
        sim2_resized = imresize(sim2, scale_factor);

        % Grayscale
        sim1_gray = rgb2gray(sim1_resized);
        sim2_gray = rgb2gray(sim2_resized);
        corners1 = detectHarrisFeatures(sim1_gray, 'MinQuality', 0.1); 
        corners2 = detectHarrisFeatures(sim2_gray, 'MinQuality', 0.1); 
        [features1, valid_points1] = extractFeatures(sim1_gray, corners1);
        [features2, valid_points2] = extractFeatures(sim2_gray, corners2);

        % Match features 
        indexPairs = matchFeatures(features1, features2);
        matchedPoints1 = valid_points1(indexPairs(:, 1), :);
        matchedPoints2 = valid_points2(indexPairs(:, 2), :);
        distances = sqrt(sum((matchedPoints1.Location - matchedPoints2.Location).^2, 2));
        movement_threshold = 10; 
        moved_indices = distances > movement_threshold;
        movedPoints1 = matchedPoints1(moved_indices);
        movedPoints2 = matchedPoints2(moved_indices);
        figure(1); clf(1);
        showMatchedFeatures(sim1_resized, sim2_resized, movedPoints1, movedPoints2, 'montage');
        

        mean_x_image_1 = mean(movedPoints1.Location(:, 1));
        mean_y_image_1 = mean(movedPoints1.Location(:, 2));
        mean_x_image_2 = mean(movedPoints2.Location(:, 1));
        mean_y_image_2 = mean(movedPoints2.Location(:, 2));
        distance_pixels = sqrt((mean_x_image_2 - mean_x_image_1)^2 + (mean_y_image_2 - mean_y_image_1)^2);

        photosize = size(sim1_gray);
        distance_pixels_ratio = distance_pixels / photosize(1, 2); 
        rgb_matrix(i, scale_index) = 100 * distance_pixels_ratio;

        elapsed_time_rgb(i, scale_index) = toc;
    end
end


figure(2); clf(2);
imagesc(rgb_matrix); 
colorbar; 
xlabel('Scale Factor Index');
ylabel('Frame Index');
xticks(1:numel(scale_factors));
xticklabels({'0.3', '0.6', '0.9', '1.2', '1.5', '1.8'});
yticks(1:17);
yticklabels(1:17);
title('Distance of object moved ');
hold on;

% Add text annotations
for i = 1:17
    for j = 1:numel(scale_factors)
        text(j, i, num2str(round(distance_matrix(i, j), 2)), ...
            'HorizontalAlignment', 'center', 'Color', 'k');
    end
end
grid on;


  %% --not used--
 FRAME_CHANGE=3;
  
 for plqual=1: size(lidar_matrix,2)
     nn=0;
    % for plfram=1: size(lidar_matrix,1)     
  for prqual=1:size(rgb_matrix,2)
    % for prfram= 1:size(rgb_matrix,1)
     
    
nn=1+nn;
     
     Differential(plqual,nn)=abs(lidar_matrix(FRAME_CHANGE,plqual) - rgb_matrix(FRAME_CHANGE,prqual) );


Spot(1,nn,plqual)=plqual;
% Spot(2,nn,plqual)=plfram;
Spot(3,nn,plqual)=prqual;
% Spot(4,nn,plqual)=prfram;


     end
 end
  %    end
  % end

 
smallest_difference= double(min(Differential, [], 2));

%[~,Pixels_equivalent(1,1)]=find(Differential==min(smallest_difference(1,1)));
[~,Pixels_equivalent(1,1)]= find(Differential==min(smallest_difference(2,1)));
[~,Pixels_equivalent(1,2)]= find(Differential==min(smallest_difference(3,1)));
[~,Pixels_equivalent(1,3)]= find(Differential==min(smallest_difference(4,1)));
[~,Pixels_equivalent(1,4)]= find(Differential==min(smallest_difference(5,1)));
[~,Pixels_equivalent(1,5)]= find(Differential==min(smallest_difference(6,1)));



%Spot_matrix(1,1)=Spot(3,Pixels_equivalent(1,1),1);
Spot_matrix(1,1)=Spot(3,Pixels_equivalent(1,1),2);
Spot_matrix(1,2)=Spot(3,Pixels_equivalent(1,2),3);
Spot_matrix(1,3)=Spot(3,Pixels_equivalent(1,3),4);
Spot_matrix(1,4)=Spot(3,Pixels_equivalent(1,4),5);
Spot_matrix(1,5)=Spot(3,Pixels_equivalent(1,5),6);

RGB_equivalence_to_lidar(1,:)=Spot(3,Spot_matrix);

lidar_range=[ 4*(8*6); 4*(16*12); 4*(32*24);4*(64*48); 4*(128*96)]; %change it to actual laserbeams
donttouch=1536*2048;
rgb_range=[1; 0.6*donttouch; 0.9*donttouch; 1.2*donttouch; 1.5*donttouch; 1.8*donttouch];

newlidar_range=lidar_range(1:5,1);
newrgb_range=rgb_range(RGB_equivalence_to_lidar(1,1:5),1);

% figure(10); clf(10);
% %hold on
% 
% loglog(newlidar_range,newrgb_range,'-r')
% xlabel('Lidar beams')
% ylabel('Megapixels (Pixels in Mega units)')




figure(11); clf(11);
hold on
plot(newlidar_range,[ 1.8874e6; 4718592;1.8874e6 ; 5.6623e6;5.6623e6],'--*r','lineWidth',1.5) %Frame 3
%plot(newlidar_range,[ 4718592; 5.6623e6; 5.6623e6; 5.6623e6;5.6623e6],'--<b','lineWidth',1.5) %Frame 12
plot(newlidar_range,[ 5.6623e6; 2.8312e6;5.6623e6; 5.6623e6;5.6623e6],'--<b','lineWidth',1.5) %Frame 14
xlim([-3000 newlidar_range(end)+2000])
xlabel('Lidar beams')
ylabel('Megapixels (Pixels in Mega units)')
legend('Frame 3','Frame 14')


%% Lidar, time response graph

figure(12); clf(12);
hold on
semilogy(newlidar_range,elapsed_time_lidar(3,2:6),'--*r','lineWidth',1.5) %Frame 3
%semilogy(newlidar_range,elapsed_time_lidar(12,2:6),'LineWidth',2) %Frame 12
semilogy(newlidar_range,elapsed_time_lidar(14,2:6),'--<b','lineWidth',1.5)  %Frame 14
xlabel('Lidar beams')
ylabel('elapsed time')
legend('period 3','period 14')


%% RGB time

rgb_range2=[ 0.6*donttouch; 0.9*donttouch; 1.2*donttouch; 1.5*donttouch; 1.8*donttouch];

figure(13); clf(13);
hold on
semilogy(rgb_range2,elapsed_time_rgb(3,2:6),'--*r','lineWidth',1.5) %Frame 3
%semilogy(rgb_range2,elapsed_time_rgb(12,2:6),'LineWidth',1) %Frame 12
semilogy(rgb_range2,elapsed_time_rgb(14,2:6),'--<b','lineWidth',1.5) %Frame 14
xlabel('Amount of Megapixels')
ylabel('elapsed time')
legend('period 3','period 14')


%% Time error - Combination of 5 


% l2-rgb2
% l3-rgb5
% l4-rgb2
% l5-rgb6
% l6-rgb6

%Frame 3
error_time1_3= abs(elapsed_time_lidar(3,2)-elapsed_time_rgb(3,2));
error_time2_3= abs(elapsed_time_lidar(3,3)-elapsed_time_rgb(3,5));
error_time3_3= abs(elapsed_time_lidar(3,4)-elapsed_time_rgb(3,2));
error_time4_3= abs(elapsed_time_lidar(3,5)-elapsed_time_rgb(3,6));
error_time5_3= abs(elapsed_time_lidar(3,6)-elapsed_time_rgb(3,6));
error_time3(:,1)=([error_time1_3;error_time2_3;error_time3_3;error_time4_3;error_time5_3]);

%Frame 14   ,[ 5.6623e6; 2.8312e6;5.6623e6; 5.6623e6;5.6623e6],
error_time1_14= abs(elapsed_time_lidar(3,2)-elapsed_time_rgb(3,6));
error_time2_14= abs(elapsed_time_lidar(3,3)-elapsed_time_rgb(3,3));
error_time3_14= abs(elapsed_time_lidar(3,4)-elapsed_time_rgb(3,6));
error_time4_14= abs(elapsed_time_lidar(3,5)-elapsed_time_rgb(3,6));
error_time5_14= abs(elapsed_time_lidar(3,6)-elapsed_time_rgb(3,6));
error_time14(:,1)=([error_time1_14;error_time2_14;error_time3_14;error_time4_14;error_time5_14]);


combinations_time=[1;2;3;4;5];

figure(14); clf(14);
hold on
plot(combinations_time,error_time3,'--*r','lineWidth',1.5)
plot(combinations_time,error_time14,'--<b','lineWidth',1.5)
xlabel('Combinations')
ylabel('Error time in seconds')
legend('Frame 3','Frame 14')


%% Distance error vs combination of 5-- not used--
%Frame 3
% error_distance1_3= abs(elapsed_time_lidar(3,2)-elapsed_time_rgb(3,2));
% error_distance2_3= abs(elapsed_time_lidar(3,3)-elapsed_time_rgb(3,5));
% error_distance3_3= abs(elapsed_time_lidar(3,4)-elapsed_time_rgb(3,2));
% error_distance4_3= abs(elapsed_time_lidar(3,5)-elapsed_time_rgb(3,6));
% error_distance5_3= abs(elapsed_time_lidar(3,6)-elapsed_time_rgb(3,6));
% error_distance3(:,1)=([error_distance1_3; error_distance2_3; error_distance3_3 ; error_distance4_3 ; error_distance5_3]);
% 
% %Frame 14   ,[ 5.6623e6; 2.8312e6;5.6623e6; 5.6623e6;5.6623e6],
% error_distance1_14= abs(elapsed_time_lidar(3,2)-elapsed_time_rgb(3,6));
% error_distance2_14= abs(elapsed_time_lidar(3,3)-elapsed_time_rgb(3,3));
% error_distance3_14= abs(elapsed_time_lidar(3,4)-elapsed_time_rgb(3,6));
% error_distance4_14= abs(elapsed_time_lidar(3,5)-elapsed_time_rgb(3,6));
% error_distance5_14= abs(elapsed_time_lidar(3,6)-elapsed_time_rgb(3,6));
% error_time14(:,1)=([error_distance1_14;error_distance2_14;error_distance3_14; error_distance4_14; error_distance5_14]);


figure(15); clf(15);
hold on
plot(combinations_time,[9.1117; 7.4675; 22.2527;0.1548; 0.5745],'--*r','lineWidth',1.5)
plot(combinations_time,[1.4118;0.0237;0.4699;0.1118;0.0334],'--<b','lineWidth',1.5)
xlabel('Combinations')
ylabel('Distance error in cm')
legend('Frame 3','Frame 14')

%% Efficiency (time_error/distance_error) over combination ---not used---

efficiency_frame3= 1/ ( (error_time3/max(error_time3)) + ([9.1117; 7.4675; 22.2527;0.1548; 0.5745]./22.2527) );
efficiency_frame14= 1/ (error_time14 + [1.4118;0.0237;0.4699;0.1118;0.0334]) ;
 % 
 % efficiency_frame33= [9.1117; 7.4675; 22.2527;0.1548; 0.5745]./error_time3;
 % efficiency_frame144= [1.4118;0.0237;0.4699;0.1118;0.0334]./error_time14;

figure(16); clf(16);
hold on
plot(error_time3,[9.1117; 7.4675; 22.2527;0.1548; 0.5745],'*r','lineWidth',3)
plot(error_time14,[1.4118;0.0237;0.4699;0.1118;0.0334],'diamondb','lineWidth',1)
xlabel('Distance error in cm')
ylabel('Time error in seconds')
legend('Frame 3','Frame 14')


% figure(17); clf(17);
% hold on
% plot(combinations_time,efficiency_frame33,'--*r','lineWidth',1.5)
% plot(combinations_time,efficiency_frame144,'--<b','lineWidth',1.5)
% xlabel('Combinations')
% ylabel('Efficiency ( time error/ distance error) ')
% legend('Frame 3','Frame 14')


%% RGB quality graph


figure(17); clf(17);
hold on
semilogy(rgb_range2,[21;24;18;20;11] ,'--*r','lineWidth',1.5) %Frame 3 real
semilogy(rgb_range2,[0;4;1;6;1] ,'--ok','lineWidth',1.5) %Frame 3 fake
semilogy(rgb_range2,[26;18;16;29;22] ,'--<b','lineWidth',1.5) %Frame 14 real
semilogy(rgb_range2,[0;2;1;3;5],'--+g','lineWidth',1.5) %Frame 14 fake
xlabel('Amount of Megapixels')
ylabel('Ammount of moving edges detected')
legend('real edges, period 3','fake edges, period 3','real edges, period 14','fake edges, period 14')


%% Lidar, quality graph

figure(18); clf(18);
hold on
semilogy(newlidar_range,[0;0;0;0;2] ,'--*r','lineWidth',1.5) %Frame 3 real
semilogy(newlidar_range,[2;2;1;1;0] ,'--ok','lineWidth',1.5) %Frame 3 fake
semilogy(newlidar_range,[1;2;2;2;3] ,'--<b','lineWidth',1.5) %Frame 14 real
semilogy(newlidar_range,[0;0;0;0;0],'--+g','lineWidth',1.5) %Frame 14 fake
xlabel('Lidar beams')
ylabel('Ammount of moving edges detected')
legend('real edges, period 3','fake edges, period 3','real edges, period 14','fake edges, period 14')
