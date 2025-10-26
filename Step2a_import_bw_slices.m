clc;
clear;
close all;
dbstop if error;

% Prompt user to select images
[fileNames, pathName] = uigetfile(...
    {'*.jpg;*.tif;*.tiff;*.png;*.gif;*.bmp;', 'All Image Files'}, ...
    'Select stack images to be processed', ...
    'MultiSelect', 'on');

% Check for file selection
if isequal(fileNames, 0)
    error('No files selected. Exiting.');
end

% Determine the number of files selected
if iscell(fileNames)
    nFiles = numel(fileNames); % Number of files when selected as cell array
else
    nFiles = 1; % Only one file was selected
    fileNames = {fileNames}; % Convert to cell array for consistency
end

% Read the first image to determine the size
firstImagePath = fullfile(pathName, fileNames{1});
img = imread(firstImagePath);
imgSize = size(img);

% Initialize a 3D matrix for the stack of images
if isa(img, 'double')
    BW3D = zeros(imgSize(1), imgSize(2), nFiles);
elseif isa(img, 'single')
    BW3D = zeros(imgSize(1), imgSize(2), nFiles, 'single');
elseif isa(img, 'uint16')
    BW3D = zeros(imgSize(1), imgSize(2), nFiles, 'uint16');
elseif isa(img, 'uint8')
    BW3D = zeros(imgSize(1), imgSize(2), nFiles, 'uint8');
else
    BW3D = false(imgSize(1), imgSize(2), nFiles);
end

% Loop through each selected file and read images into the 3D matrix
for iFile = 1:nFiles
    currentFilePath = fullfile(pathName, fileNames{iFile}); % Get the current file path
    if isfile(currentFilePath) % Ensure the file exists
        BW3D(:,:,iFile) = imread(currentFilePath);
    else
        error(['File not found: ' currentFilePath]);
    end
end

% Optional: Binarize the 3D matrix slice by slice
BW = false(size(BW3D)); % Initialize binarized matrix

% Use double conversion for binarizing
for iSlice = 1:nFiles
    % Ensure the data type is suitable for binarization
    if islogical(BW3D(:,:,iSlice))
        % If it's already logical, just assign it
        BW(:,:,iSlice) = BW3D(:,:,iSlice); 
    else
        % Convert to double, if necessary and then binarize
        BW(:,:,iSlice) = imbinarize(double(BW3D(:,:,iSlice))); % Binarize each slice
    end
end

% Get the first object's binary image from region properties
stats = regionprops3(BW, 'Image', 'BoundingBox');
if ~isempty(stats)
    BW = stats.Image{1, 1}; % Extract the first object's binary image
else
    error('No objects found in the binary volume.');
end

% Save the processed binary matrix
save('0_90_40kV_250uA_16um_res_rec.mat', 'BW');

% Visualization of the 3D binary matrix
imshow3Dfull(BW);