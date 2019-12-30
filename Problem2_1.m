clear; clc;

Problem_2_1();

function Problem_2_1()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
Xdir = uigetdir('D:\KAIST\Courses\dip\hw\hw4\Test_images');
file = fopen(fullfile(Xdir,'\Gray_monarch_512x512.raw'),'rb');
gray_monarch = fread(file,fliplr([512,512]),'*uint8')';
fclose(file);

file = fopen(fullfile(Xdir,'\Color_monarch_512x512.raw'),'rb');
color_monarch = fread(file,fliplr([512,512*3]),'*uint8')';
fclose(file);
%%
%%---------------------Insert code below ----------------------%%
gray_monarch = double(gray_monarch);
color_monarch = double(color_monarch);

% Problem 2.1a - Histogram Equalization
gray_hist = histEq(gray_monarch);

[~,~,~,rgb_monarch] = divRGB(color_monarch);
rgb_hist = zeros(512,512,3);
rgb_hist(:,:,1) = histEq(rgb_monarch(:,:,1));
rgb_hist(:,:,2) = histEq(rgb_monarch(:,:,2));
rgb_hist(:,:,3) = histEq(rgb_monarch(:,:,3));

[~,~,~,hsi_monarch] = RGB2HSI(rgb_monarch);
hsi_monarch = uint8(hsi_monarch*255);
hsi_hist = zeros(512,512,3);
hsi_hist(:,:,1) = histEq(hsi_monarch(:,:,1))/255;
hsi_hist(:,:,2) = histEq(hsi_monarch(:,:,2))/255;
hsi_hist(:,:,3) = histEq(hsi_monarch(:,:,3))/255;

[~,~,~,cmy_monarch] = RGB2CMY(rgb_monarch);
cmy_monarch = uint8(cmy_monarch);
cmy_hist = zeros(512,512,3);
cmy_hist(:,:,1) = histEq(cmy_monarch(:,:,1));
cmy_hist(:,:,2) = histEq(cmy_monarch(:,:,2));
cmy_hist(:,:,3) = histEq(cmy_monarch(:,:,3));

% Problem 2.1b - Histogram Modification
alpha1 = 0.3; alpha2 = 0.8;
gray_alpha1 = histEqMod(gray_monarch,alpha1);
gray_alpha2 = histEqMod(gray_monarch,alpha2);

%% Displaying figures (Edit this part as needed)
figure('Name', 'Problem 2.1a - Histogram Equalization');
subplot(2,4,1); imshow(uint8(gray_monarch)); title('Gray Monarch');
subplot(2,4,2); imshow(uint8(gray_hist)); title('Gray Histogram Equalization');
subplot(2,4,5); imshow(uint8(rgb_monarch)); title('Color Monarch');
subplot(2,4,6); imshow(uint8(rgb_hist)); title('RGB Histogram Equalization');
subplot(2,4,7); imshow(uint8(255*hsi_hist)); title('HSI Histogram Equalization');
subplot(2,4,8); imshow(uint8(cmy_hist)); title('CMY Histogram Equalization');

figure('Name', 'Problem 2.1b - Histogram Modification');
subplot(1,3,1); imshow(uint8(gray_monarch)); title('Gray Monarch');
subplot(1,3,2); imshow(uint8(gray_alpha1)); title('alpha = 0.3');
subplot(1,3,3); imshow(uint8(gray_alpha2)); title('alpha = 0.8');

%%---------------------------------------------------------------%%
end

%% Histogram Equalization
function A = histEq(X)
	[m, n] = size(X);

	% Get pdf
	P = zeros(256,1);
	for i = 1:256
		P(i) = sum(sum(X==i-1))/(m*n);
	end

	% Get cdf
	C = zeros(256,1);
	C(1) = P(1);
	for i = 2:256
		C(i) = C(i-1) + P(i);
	end

	% Mapping
	T = C * 255;
	A = zeros(m,n);
	for i = 1:m
		for j = 1:n
			A(i,j) = T(uint8(X(i,j))+1);
		end
	end
end

%% Histogram Equalization Modification
function A = histEqMod(X,alpha)
	[m, n] = size(X);

	% Get pdf
	P = zeros(256,1);
	for i = 1:256
		P(i) = sum(sum(X==i-1))/(m*n);
	end

	% Get cdf
	C = zeros(256,1);
	C(1) = P(1);
	for i = 2:256
		C(i) = C(i-1) + P(i);
	end

	% Mapping
	A = -log(1-C(X+1))/alpha;
	A = round(A*255);
end

%% Separating RGB components
function [rA, gA, bA, A] = divRGB(X)
	[m, n] = size(X);
	p = n / 3;

	rA = zeros(m,p,3); rA(:,:,1) = X(:,1:3:(n-2)); rA = uint8(rA);  % Red
	gA = zeros(m,p,3); gA(:,:,2) = X(:,2:3:(n-1)); gA = uint8(gA);  % Green
	bA = zeros(m,p,3); bA(:,:,3) = X(:,3:3:n); bA = uint8(bA);      % Blue

	A = rA + gA + bA;    % RGB image
end

%% Converting RGB to HSI
function [hA, sA, iA, A] = RGB2HSI(X)
	[m, n, p] = size(X);
	A = zeros(m, n, p);

	X_norm = double(X)/255;			% Normalize input RGB image

	rX_norm = X_norm(:,:,1);    % Red normalized
	gX_norm = X_norm(:,:,2);    % Green normalized
	bX_norm = X_norm(:,:,3);    % Blue normalized

	iA = (rX_norm + gX_norm + bX_norm)/3;   % Intensity
	sA = 1 - min(X_norm,[],3)./iA;          % Saturation

	% Calculate hue
	V1 = 2 * rX_norm - gX_norm - bX_norm;
	V2 = 2 * sqrt((rX_norm - gX_norm).^2 + (rX_norm - bX_norm).*(gX_norm - bX_norm));
	theta = acos(V1./V2);
	hA = theta/(2*pi);

	A(:,:,1) = hA;
	A(:,:,2) = sA;
	A(:,:,3) = iA;
end

%% Converting RGB to CMY
function [cA, mA, yA, A] = RGB2CMY(X)
	[m, n, p] = size(X);
	A = zeros(m, n, p);

	X_norm = double(X)/255;			% Normalize input RGB image

	rX_norm = X_norm(:,:,1);    % Red normalized
	gX_norm = X_norm(:,:,2);    % Green normalized
	bX_norm = X_norm(:,:,3);    % Blue normalized

	cA = 1 - rX_norm;
	mA = 1 - gX_norm;
	yA = 1 - bX_norm;

	A(:,:,1) = cA * 255;
	A(:,:,2) = mA * 255;
	A(:,:,3) = yA * 255;
end
