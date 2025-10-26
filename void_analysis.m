%% Local Porosity and Void Analysis Calculation
clear; clc; close all;

% Load the 3D binary matrix
type = '0_90_40kV_250uA_16um_res_rec'; % Change to the desired type
load([type, '.mat']); % Assumed BW is loaded here

% Parameters
n = 5; % Kernel size (must be an odd number)
window_center = (n + 1) / 2;

% Prepare scanning window matrix and output matrix
sz = size(BW);
output_porosity = zeros(sz(1) - n + 1, sz(2) - n + 1, sz(3) - n + 1);
output_void_diameter = zeros(size(output_porosity));

%% Traverse the entire 3D binary matrix to calculate local porosity and voids
for i = window_center:sz(1)-window_center
    for j = window_center:sz(2)-window_center
        for k = window_center:sz(3)-window_center

            % Extracting the local scanning window
            scan_window = BW(i-(window_center-1):i+(window_center-1), ...
                              j-(window_center-1):j+(window_center-1), ...
                              k-(window_center-1):k+(window_center-1));
            
            % Local porosity calculation
            local_porosity = 1 - (sum(scan_window, 'all') / n^3);
            output_porosity(i-(window_center-1), j-(window_center-1), k-(window_center-1)) = local_porosity;
            
            % Example of void size calculation (diameter based on scanned window)
            % Here, we're assuming the voids are the "0"s in BW and material "1"s
            void_voxels = sum(scan_window(:) == 0); % Count zeros
            if void_voxels > 0  % To avoid division by zero
                % Calculate an approximate diameter; this can be refined
                void_diameter = (6 * void_voxels / (4/3 * pi))^(1/3);  % Based on volume of a sphere
                output_void_diameter(i-(window_center-1), j-(window_center-1), k-(window_center-1)) = void_diameter;
            end
        end
    end
end

%% Save the results
save([type, '_local_porosity_', num2str(n), '.mat'], 'output_porosity', 'output_void_diameter', '-v7.3');

%% Post-Processing Results, Example Visualization
% Visualize a slice of the porosity
% figure;
% slice = round(size(output_porosity, 3) / 2); % Middle slice
% imagesc(output_porosity(:,:,slice));
% colorbar;
% title('Local Porosity Slice');
% xlabel('X-axis');
% ylabel('Y-axis');
% axis image;

figure;
slice = round(size(output_void_diameter, 3) / 2);  % Choose middle slice
imagesc(output_void_diameter(:,:,slice));
axis image;
colorbar;
caxis([0, max(output_void_diameter(:))]); % Optional: fix color scale
title('Local Void Diameter (Voxel Units)');
xlabel('X-axis');
ylabel('Y-axis');

cb = colorbar;
cb.Label.String = 'Diameter (voxels)';

% Flatten and remove zeros
diameters = output_void_diameter(:);
diameters = diameters(diameters > 0);

% Histogram
figure;
histogram(diameters, 30);  % Adjust bin count as needed
xlabel('Void Diameter (voxels)');
ylabel('Frequency');
title('Histogram of Local Void Diameters');



% Example of analyzing the orientation of voids (for a specific slice)
% Orientation analysis could be added here based on specific criteria