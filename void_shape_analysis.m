function void_shape_analysis(BW)
% void_shape_analysis - Analyze sphericity and compactness of 3D voids in BW
%   void_shape_analysis(BW) analyzes each void in the 3D binary BW (1=material, 0=void)
%   void_shape_analysis      (no argument) looks for BW in workspace.
%
% Computes per-void shape descriptors: sphericity and compactness.
% Sphericity is bounded in [0,1]; compactness is a normalized volume-to-length ratio.

% --- Config ---
voxelSize = [1, 1, 1]; % Set to [dx dy dz] if anisotropic voxels
showNExamples = 6;

% --- Input Handling ---
if nargin < 1
    if evalin('base','exist(''BW'',''var'')')
        BW = evalin('base','BW');
        disp('Using BW variable from workspace');
    else
        error('No BW input given, and no variable "BW" found in workspace.');
    end
end

if ~islogical(BW)
    BW = logical(BW);
end
if ndims(BW) ~= 3
    error('BW must be a 3D logical (binary) array.');
end

voids = ~BW; % voids are 1, material is 0
CC = bwconncomp(voids, 26);
if CC.NumObjects == 0
    error('No voids detected in BW!');
end

% --- Measure Properties ---
props = regionprops3(voids, 'Volume','SurfaceArea','PrincipalAxisLength','VoxelList','Centroid');

% Scale volumes/surface area if voxel size is not isotropic
V = props.Volume .* prod(voxelSize);       % Volume in real units
S = props.SurfaceArea .* (mean(voxelSize)^2); % Approximate surface area in real units

% --- Compute Shape Descriptors ---
% Remove invalid entries
valid = S > 0 & V > 0;
V = V(valid);
S = S(valid);
principalAxes = props.PrincipalAxisLength(valid,:);
voxelLists = props.VoxelList(valid);

% Sphericity (bounded in [0,1])
sphericity = (pi^(1/3)) * (6 * V).^(2/3) ./ S;
sphericity(sphericity > 1 | sphericity < 0) = NaN;

% Compactness (normalized volume-to-length ratio)
max_axis = max(principalAxes, [], 2)/2 + eps; % Radius approximation
compactness = (V.^(1/3)) ./ max_axis;
compactness(compactness < 0 | ~isfinite(compactness)) = NaN;

% --- Plot Distributions ---
figure('Name','Void Shape Analysis');
subplot(1,2,1)
histogram(sphericity(~isnan(sphericity)), 30, 'FaceColor', [0 0.6 1]);
xlabel('Sphericity'); ylabel('Count');
title('Void Sphericity Distribution');

subplot(1,2,2)
histogram(compactness(~isnan(compactness)), 30, 'FaceColor', [1 0.4 0.1]);
xlabel('Compactness'); ylabel('Count');
title('Void Compactness Distribution');

sgtitle('Void Shape Descriptors');

% --- Show Example Voids (colored by sphericity) ---
N = min(numel(V), showNExamples);
[~, sortIdx] = sort(sphericity, 'descend');
exIdxs = sortIdx(1:N);

figure('Name','Example Voids: 3D Scatter');
for i = 1:N
    subplot(2, ceil(N/2), i);
    this_vox = voxelLists{exIdxs(i)} .* voxelSize; % scale to real units
    scatter3(this_vox(:,1), this_vox(:,2), this_vox(:,3), 10, ...
        sphericity(exIdxs(i))*ones(size(this_vox,1),1), '.');
    axis equal tight
    colormap parula
    colorbar
    title(sprintf('Void %d:\nSph = %.2f', exIdxs(i), sphericity(exIdxs(i))));
    xlabel('X'); ylabel('Y'); zlabel('Z');
end
sgtitle('Example Voids (Each Scatter3 Colored by Sphericity)');

end



% function void_shape_analysis(BW)
% % void_shape_analysis - Analyze sphericity and compactness of 3D voids in BW
% %   void_shape_analysis(BW) analyzes each void in the 3D binary BW (1=material, 0=void)
% %   void_shape_analysis      (no argument) looks for BW in workspace.
% %
% % This function computes sphericity and compactness for each void, and
% % visualizes their distributions and shows example 3D voids.
% 
% % --- Robust Input Handling ---
% if nargin < 1
%     if evalin('base','exist(''BW'',''var'')')
%         BW = evalin('base','BW');
%         disp('Using BW variable from workspace');
%     else
%         error('No BW input given, and no variable "BW" found in workspace.');
%     end
% end
% 
% if ~islogical(BW)
%     BW = logical(BW);
% end
% if ndims(BW) ~= 3
%     error('BW must be a 3D logical (binary) array.');
% end
% 
% voids = ~BW; % now 1=void, 0=material
% CC = bwconncomp(voids, 26);
% if CC.NumObjects == 0
%     error('No voids detected in BW!');
% end
% 
% props = regionprops3(voids, 'Volume','SurfaceArea','PrincipalAxisLength','VoxelList','Centroid');
% 
% V = props.Volume;      % in voxels
% S = props.SurfaceArea; % in voxels
% sphericity = (pi^(1/3))*(6.*V).^(2/3)./S; % Standard formula
% 
% max_axis = max(props.PrincipalAxisLength,[],2)/2 + eps; % +eps avoids zero-division
% compactness = V ./ (max_axis.^3);
% 
% figure('Name','Void Shape Analysis');
% subplot(1,2,1)
% histogram(sphericity,30,'FaceColor',[0 0.6 1]);
% xlabel('Sphericity'); ylabel('Count');
% title('Void Sphericity Distribution');
% subplot(1,2,2)
% histogram(compactness,30,'FaceColor',[1 0.4 0.1]);
% xlabel('Compactness'); ylabel('Count');
% title('Void Compactness Distribution');
% sgtitle('Void Shape Descriptors');
% 
% % --- Show some example voids in 3D, colored by sphericity ---
% N = min(numel(V), 6); % Show up to 6
% [~,sortIdx] = sort(sphericity,'descend');
% exIdxs = sortIdx(1:N);
% 
% figure('Name','Example Voids: 3D scatter');
% for i = 1:N
%     subplot(2, ceil(N/2), i);
%     this_vox = props.VoxelList{exIdxs(i)};
%     scatter3(this_vox(:,1), this_vox(:,2), this_vox(:,3), 10, ...
%         sphericity(exIdxs(i))*ones(size(this_vox,1),1), '.');
%     axis equal tight
%     colormap parula
%     colorbar
%     title(sprintf('Void %d:\nSph=%.2f',exIdxs(i),sphericity(exIdxs(i))));
%     xlabel('X'); ylabel('Y'); zlabel('Z');
% end
% sgtitle('Example Voids (Each Scatter3 colored by Sphericity)');