clear; clc;

Problem_2_3();

function Problem_2_3()
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
Gau_noise = gaussianNoise(gray_monarch,0,0.01);

%% DFT
Ldft = zonal_mask(gray_monarch,80,'L','DFT');
Hdft = zonal_mask(gray_monarch,80,'H','DFT');

% No noise
fft2D_monarch = FFT(FFT(gray_monarch).').';
recLdft = IFFT(IFFT(Ldft.*fft2D_monarch).').';
recHdft = IFFT(IFFT(Hdft.*fft2D_monarch).').';

% Gaussian noise
fft2D_Gau = FFT(FFT(Gau_noise).').';
recLdft_Gau = IFFT(IFFT(Ldft.*fft2D_Gau).').';
recHdft_Gau = IFFT(IFFT(Hdft.*fft2D_Gau).').';

%% DCT
Ldct = zonal_mask(gray_monarch,80,'L','DCT');
Hdct = zonal_mask(gray_monarch,80,'H','DCT');

% No noise
dct2D_monarch = DCT2D(gray_monarch,1);
recLdct = IDCT2D(Ldct.*dct2D_monarch,1);
recHdct = IDCT2D(Hdct.*dct2D_monarch,1);

% Gaussian noise
dct2D_Gau = DCT2D(Gau_noise,1);
recLdct_Gau = IDCT2D(Ldct.*dct2D_Gau,1);
recHdct_Gau = IDCT2D(Hdct.*dct2D_Gau,1);

%% Displaying figures (Edit this part as needed)
figure('Name','Problem 2.3 - Transform Domain Filtering');
subplot(2,4,1); imshow(uint8(gray_monarch)); title('Gray Monarch');
subplot(2,4,2); imshow(uint8(Gau_noise)); title('Gaussian Noise');
subplot(2,4,3); imshow(uint8(abs(recLdft))); title('LPF no noise DFT');
subplot(2,4,4); imshow(uint8(abs(recHdft))); title('HPF no noise DFT');
subplot(2,4,5); imshow(uint8(abs(recLdft_Gau))); title('LPF Gaussian DFT');
subplot(2,4,6); imshow(uint8(abs(recHdft_Gau))); title('HPF Gaussian DFT');
subplot(2,4,7); imshow(uint8(abs(recLdct_Gau))); title('LPF Gaussian DCT');
subplot(2,4,8); imshow(uint8(abs(recHdct_Gau))); title('HPF Gaussian DCT');

%%---------------------------------------------------------------%%
end

%% FFT Function
function A = FFT(X)
  N = length(X);

  X_even = X(:,1:2:N-1);
  X_odd = X(:,2:2:N);
  W = exp(-1j*4*pi/N);

  m = 0:1:(N/2-1);
  k = 0:1:(N-1);

  KM_even = m' * k;
  KM_odd = (m+0.5)' * k;

  WN_even = W.^KM_even;
  WN_odd = W.^KM_odd;

  A = X_even*WN_even + X_odd*WN_odd;
end

%% Inverse FFT Function
function A = IFFT(X)
  N = length(X);

  X_even = X(:,1:2:N-1);
  X_odd = X(:,2:2:N);
  W = exp(1j*4*pi/N);

  m = 0:1:(N/2-1);
  k = 0:1:(N-1);

  KM_even = m' * k;
  KM_odd = (m+0.5)' * k;

  WN_even = W.^KM_even;
  WN_odd = W.^KM_odd;

  A = (X_even*WN_even + X_odd*WN_odd)/N;
end

%% 2-Dimension DCT Function
function A = DCT2D(X,blkSize)
  n = length(X);
  if blkSize == 1     % Normal DCT
    DCT = sqrt(2/n)*cos(pi/(2*n)*(0:1:n-1)'*(1:2:2*n-1));
    DCT(1,:) = DCT(1,:)/sqrt(2);
    A = DCT * X * DCT';
  else    % Blockwise DCT
    N = n / blkSize;
    DCT = sqrt(2/N)*cos(pi/(2*N)*(0:1:N-1)'*(1:2:2*N-1));
    DCT(1,:) = DCT(1,:) / sqrt(2);
   
    BLK = zeros(N,N,blkSize,blkSize);
    for i = 0:(blkSize-1)
        for j = 0:(blkSize-1)
          BLK(:,:,i+1,j+1) = X(1+i*N:N+i*N,1+j*N:N+j*N);
          BLK(:,:,i+1,j+1) = DCT*BLK(:,:,i+1,j+1)*DCT';
        end 
    end
   
    A = size(X);
    for i = 0:(blkSize-1)
      for j = 0:(blkSize-1)
        A(1+i*N:N+i*N,1+j*N:N+j*N) = BLK(:,:,i+1,j+1);
      end
    end
  end
end

%% Inverse 2-Dimension DCT Function (IDCT)
function A = IDCT2D(X,blkSize)
  n = length(X);
  if blkSize == 1     % Normal IDCT
    IDCT = sqrt(2/n)*cos(pi/(2*n)*(1:2:2*n-1)'*(0:1:n-1));
    IDCT(1,:) = IDCT(1,:)/sqrt(2);
    A = IDCT * X * IDCT';
  else
    N = n / blkSize;
    IDCT = sqrt(2/N)*cos(pi/(2*N)*(1:2:2*N-1)'*(0:1:N-1));
    IDCT(1,:) = IDCT(1,:) / sqrt(2);
   
    BLK = zeros(N,N,blkSize,blkSize);
    for i = 0:(blkSize-1)
        for j = 0:(blkSize-1)
          BLK(:,:,i+1,j+1) = X(1+i*N:N+i*N,1+j*N:N+j*N);
          BLK(:,:,i+1,j+1) = IDCT * BLK(:,:,i+1,j+1) * IDCT';
        end
    end
   
    A = size(X);
    for i = 0:(blkSize-1)
        for j = 0:(blkSize-1)
          A(1+i*N:N+i*N,1+j*N:N+j*N) = BLK(:,:,i+1,j+1); 
        end
    end
  end
end

function A = gaussianNoise(X,mu,sigma)
	A = X/255 + sqrt(sigma)*randn(size(X)) + mu;
	A = A * 255;
end

function A = zonal_mask(X,m,filter_type,transform_type) 
	% n : whole size of mask (same as image)
	% m : range of filter
	% filter_type: L (LPF) or H (HPF)
	% transform_type: DCT or DFT
	[M, N] = size(X);

  if transform_type == 'DFT'
    n = M;
    A = zeros(n,n);
    A_tmp = zeros(n,n);
    A_tmp1 = zeros(n,n);
    A_tmp2 = zeros(n,n);
    A_tmp3 = zeros(n,n);

    for i = 1:m
      for j = 1:m
        A_tmp(i,j) = (j < round(sqrt(m^2-i^2)));  % For circle shape
      end
    end

    for i = 1:m
      for j = 1:m
        A_tmp1(n-i+1,n-j+1) = A_tmp(i,j);
        A_tmp2(n+1-i,j) = A_tmp(i,j);
        A_tmp3(i,n-j+1) = A_tmp(i,j);
      end
    end
    A = A_tmp + A_tmp1 + A_tmp2 + A_tmp3;

    if filter_type == 'H'
      A = ones(n,n) - A;
    end
  end

  if transform_type == 'DCT'
    n = M;
    A = zeros(n,n);
    A_tmp = zeros(n,n);

    for i = 1:m
      for j = 1:m
        A_tmp(i,j) = (j < round(sqrt(m^2-i^2)));  % For circle shape
      end
    end

    A = A_tmp;

    if filter_type == 'H'
      A = ones(n,n) - A;
    end
  end
end
