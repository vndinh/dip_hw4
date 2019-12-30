clear; clc;

Problem_2_2();

function Problem_2_2()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
Xdir = uigetdir('D:\KAIST\Courses\dip\hw\hw4\Test_images');
file = fopen(fullfile(Xdir,'\Gray_monarch_512x512.raw'),'rb');
gray_monarch = fread(file,fliplr([512,512]),'*uint8')';
fclose(file);
%%
%%---------------------Insert code below ----------------------%%
gray_monarch = double(gray_monarch);

% Problem 2.2a - Gamma Correction
gamma1 = 0.45; gamma2 = 2.2;
gray_gamma1 = gammaCorrection(gray_monarch,gamma1);
gray_gamma2 = gammaCorrection(gray_monarch,gamma2);

% Problem 2.2b - Spartial Filtering
SP_noise = saltPepperNoise(gray_monarch,4,251);
Gau_noise = gaussianNoise(gray_monarch,0,0.01);

SP_3_mean = MeanFilter(SP_noise,3);
SP_7_mean = MeanFilter(SP_noise,7);
SP_3_med = MedianFilter(SP_noise,3);
SP_7_med = MedianFilter(SP_noise,7);
SP_3_dir = DirSmooth(SP_noise,3);
SP_7_dir = DirSmooth(SP_noise,7);

Gau_3_mean = MeanFilter(Gau_noise,3);
Gau_7_mean = MeanFilter(Gau_noise,7);
Gau_3_med = MedianFilter(Gau_noise,3);
Gau_7_med = MedianFilter(Gau_noise,7);
Gau_3_dir = DirSmooth(Gau_noise,3);
Gau_7_dir = DirSmooth(Gau_noise,7);

PSNR_SP_3_mean = PSNR(gray_monarch,SP_3_mean);
PSNR_SP_7_mean = PSNR(gray_monarch,SP_7_mean);
PSNR_SP_3_med = PSNR(gray_monarch,SP_3_med);
PSNR_SP_7_med = PSNR(gray_monarch,SP_7_med);
PSNR_SP_3_dir = PSNR(gray_monarch,SP_3_dir);
PSNR_SP_7_dir = PSNR(gray_monarch,SP_7_dir);

PSNR_Gau_3_mean = PSNR(gray_monarch,Gau_3_mean);
PSNR_Gau_7_mean = PSNR(gray_monarch,Gau_7_mean);
PSNR_Gau_3_med = PSNR(gray_monarch,Gau_3_med);
PSNR_Gau_7_med = PSNR(gray_monarch,Gau_7_med);
PSNR_Gau_3_dir = PSNR(gray_monarch,Gau_3_dir);
PSNR_Gau_7_dir = PSNR(gray_monarch,Gau_7_dir);

%% Displaying figures (Edit this part as needed)
figure('Name', 'Problem 2.2a - Gamma Correction');
subplot(1,3,1); imshow(uint8(gray_monarch)); title('Gray Monarch');
subplot(1,3,2); imshow(uint8(gray_gamma1)); title('gamma = 0.45');
subplot(1,3,3); imshow(uint8(gray_gamma2)); title('gamma = 2.2');

figure('Name', 'Problem 2.2b - Spatial Filtering');
subplot(4,4,1); imshow(uint8(gray_monarch)); title('Gray Monarch');
subplot(4,4,2); imshow(uint8(SP_noise)); title('Salt & Pepper Noise');
subplot(4,4,3); imshow(uint8(Gau_noise)); title('Gaussian Noise');
subplot(4,4,4); imshow(uint8(SP_3_mean)); title('SP 3 mean');
subplot(4,4,5); imshow(uint8(Gau_3_mean)); title('Gau 3 mean');
subplot(4,4,6); imshow(uint8(SP_7_mean)); title('SP 7 mean');
subplot(4,4,7); imshow(uint8(Gau_7_mean)); title('Gau 7 mean');
subplot(4,4,8); imshow(uint8(SP_3_med)); title('SP 3 med');
subplot(4,4,9); imshow(uint8(Gau_3_med)); title('Gau 3 med');
subplot(4,4,10); imshow(uint8(SP_7_med)); title('SP 7 med');
subplot(4,4,11); imshow(uint8(Gau_7_med)); title('Gau 7 med');
subplot(4,4,12); imshow(uint8(SP_3_dir)); title('SP 3 dir');
subplot(4,4,13); imshow(uint8(Gau_3_dir)); title('Gau 3 dir');
subplot(4,4,14); imshow(uint8(SP_7_dir)); title('SP 7 dir');
subplot(4,4,15); imshow(uint8(Gau_7_dir)); title('Gau 7 dir');

fprintf('PSNR_SP_3_mean = %f \n', PSNR_SP_3_mean);
fprintf('PSNR_SP_7_mean = %f \n', PSNR_SP_7_mean);
fprintf('PSNR_SP_3_med = %f \n', PSNR_SP_3_med);
fprintf('PSNR_SP_7_med = %f \n', PSNR_SP_7_med);
fprintf('PSNR_SP_3_dir = %f \n', PSNR_SP_3_dir);
fprintf('PSNR_SP_7_dir = %f \n', PSNR_SP_7_dir);

fprintf('PSNR_Gau_3_mean = %f \n', PSNR_Gau_3_mean);
fprintf('PSNR_Gau_7_mean = %f \n', PSNR_Gau_7_mean);
fprintf('PSNR_Gau_3_med = %f \n', PSNR_Gau_3_med);
fprintf('PSNR_Gau_7_med = %f \n', PSNR_Gau_7_med);
fprintf('PSNR_Gau_3_dir = %f \n', PSNR_Gau_3_dir);
fprintf('PSNR_Gau_7_dir = %f \n', PSNR_Gau_7_dir);

%%---------------------------------------------------------------%%
end

%% Gamma Correction
function A = gammaCorrection(X,gamma)
	A = (X/255).^gamma*255;
end

function A = saltPepperNoise(X,black,white)
	[m, n] = size(X);
	A = X;
	B = randi(255,m,n);
	A(B <= black) = 0;
	A(B >= white) = 255;
end

function A = gaussianNoise(X,mu,sigma)
	A = X/255 + sqrt(sigma)*randn(size(X)) + mu;
	A = A * 255;
end

function A = conv2Dsame(X,H)
	[xRows, xCols] = size(X);
	[hRows, hCols] = size(H);

	% Rotate kernel matrix by 180 degree
	K = rot90(H,2);

	% Zero padding matrix
	center = floor((size(K)+1)/2);
	left = center(2) - 1;
	right = hCols - center(2);
	top = center(1) - 1;
	bot = hRows - center(1);
	Z = zeros(xRows + top + bot, xCols + left + right);

	for i = (1+top) : (xRows+top)
    for j = (1+left) : (xCols+left)
      Z(i,j) = X(i - top, j - left);
    end
	end

	% Convolution
	A = zeros(xRows,xCols);
	for m = 1:xRows
   	for n = 1:xCols
      for i = 1:hRows
        for j = 1:hCols
          p = m - 1;
          q = n - 1;
          A(m,n) = A(m,n) + Z(i + p,j + q) * K(i,j);
        end
      end
   	end
	end
end

function A = MeanFilter(X,N)
	H = ones(N,N)/N^2;
	A = conv2Dsame(X,H);
end

function A = MedianFilter(X,N)
	[m, n] = size(X);
	winSize = N*N;

	% Zero padding
	p = floor(N/2);
	Padd = zeros(m+N-1,n+N-1);
	for i = 1:m
		for j = 1:n
			Padd(i+p,j+p) = X(i,j);
		end
	end

	W = zeros(winSize);
	for x = (p+1):(m+p)
		for y = (p+1):(n+p)
			i = 1;
			for wx = 0:N-1
				for wy = 0:N-1
					W(i) = Padd(x+wx-p,y+wy-p);
					i = i + 1;
				end
			end
			W = sort(W);
			A(x-p,y-p) = W(floor(winSize/2));
		end
	end
end

function A = DirSmooth(X,N)
  [m, n] = size(X);
  p = floor(N/2); % Kernel range (positive range number)

  Padd = zeros(m+N-1,n+N-1);  % Zero Padding
  for i = 1:m
    for j = 1:n          
      Padd(i+p,j+p) = X(i,j);            
    end
  end

  A = ones(m,n);
  for i = 1:m
    for j = 1:n
      temp = zeros(1,N); % for window
      degree_vec = zeros(1,4); % for 0, 45, 90, 135 degree

      % 0 degree
      % count number to add pixel value in window into temp vector
      count = 1;
      for t = -p:p
        temp(count) = Padd(i+p,j+p+t);
        count = count + 1;
      end
      degree_vec(1)=sum(temp)/N; % average for pixel value in window

      % 45 degree
      % Count number to add pixel value in window into temp vector
      count = 1;
      for t = -p:p
      	temp(count) = Padd(i+p-t,j+p+t);
        count = count + 1;
      end
      degree_vec(2)=sum(temp)/N; % Average value in window

      % 90 degree
      % Count number to add pixel value in window into temp vector
      count = 1;
      for t = -p:p
        temp(count)=Padd(i+p+t,j+p);
        count = count + 1;
      end
      degree_vec(3) = sum(temp)/N; % Average value in window

      % 135 degree
      % Count number to add pixel value in window into temp vector
      count = 1;
      for t = -p:p
        temp(count) = Padd(i+p+t,j+p+t);
        count = count + 1;
      end
      degree_vec(4) = sum(temp)/N; % Average value in window
            
      compare_vec = Padd(i+p,j+p)*ones(1,4);
      select_vec = abs(compare_vec-degree_vec);
      [~,idx] = min(select_vec);
            
      A(i,j) = degree_vec(idx); % Choose the closest direction value
    end
  end
end

function psnr = PSNR(Orig,Dist)
	[m, n, p] = size(Orig);
	Orig = double(Orig);
	Dist = double(Dist);
	error = Orig - Dist;
	MSE = sum(sum(sum(error.^2)))/(m*n*p);
	if MSE > 0
    psnr = 20*log10(max(max(max(Orig)))) - 10*log10(MSE);
	else
    psnr = 99;
	end
end
