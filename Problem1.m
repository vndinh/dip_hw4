clear; clc;

Problem_1();

function Problem_1()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
imgdir = uigetdir('D:\KAIST\Courses\dip\hw\hw4\Test_images');
file = fopen(fullfile(imgdir,'\Gray_monarch_512x512.raw'),'rb');
gray_monarch = fread(file,fliplr([512,512]),'*uint8')';
fclose(file);
%%
%%---------------------Insert code below ----------------------%%
gray_monarch = double(gray_monarch);
%% Problem 1.1 - DFT
% FFT
fft2D_monarch = FFT(FFT(gray_monarch).').';
Re_fft2D = real(fft2D_monarch);
Im_fft2D = imag(fft2D_monarch);
fft_mag = abs(fft2D_monarch);
log_fft_mag = 10*log(fft_mag);
fft_phase = atan2(Im_fft2D,Re_fft2D);
disp_fft_phase = (fft_phase+pi)*256/(2*pi);

% Inverse FFT
ifft2D_monarch = IFFT(IFFT(fft2D_monarch).').';
ifft2D_magnitude = IFFT(IFFT(fft_mag).').';
fft_phase_complex = 256*exp(1j*fft_phase);
ifft2D_phase = IFFT(IFFT(fft_phase_complex).').';

disp_ifft_mag = uint8(abs(ifft2D_magnitude));
disp_ifft_phase = uint8(abs(ifft2D_phase)*512);

%% Problem 1.2 - DCT
% DCT
dct_monarch = DCT2D(gray_monarch,1);
dct16_monarch = DCT2D(gray_monarch,16);
dct_mag = 10*log(abs(dct_monarch));
dct16_mag = 10*log(abs(dct16_monarch));

% IDCT
idct_monarch = IDCT2D(dct_monarch,1);
idct16_monarch = IDCT2D(dct16_monarch,16);

%% Problem 1.3 - DHT
dht_monarch = DHT(gray_monarch);
idht_monarch = IDHT(dht_monarch);
dht_mag = 30*log(abs(dht_monarch));

%% Problem 1.5 - DWT
% 3-level DWT
dwt_1 = DWT2D(DWT2D(gray_monarch).').';
dwt_2 = DWT2D(DWT2D(dwt_1(1:256,1:256)).').';
dwt_3 = DWT2D(DWT2D(dwt_2(1:128,1:128)).').';
dwt_1(1:256,1:256) = dwt_2;
dwt_1(1:128,1:128) = dwt_3;
dwt_1(1:64,1:64) = dwt_3(1:64,1:64);

idwt_monarch = dwt_1;
idwt_monarch(1:64,1:64) = idwt_monarch(1:64,1:64);
idwt_monarch(1:128,1:128) = IDWT2D(IDWT2D(idwt_monarch(1:128,1:128)).').';
idwt_monarch(1:256,1:256) = IDWT2D(IDWT2D(idwt_monarch(1:256,1:256)).').';
idwt_monarch = IDWT2D(IDWT2D(idwt_monarch).').';

%% Displaying figures (Edit this part as needed)
% Problem 1.1
figure('Name', 'Problem 1.1 - DFT');
subplot(2,3,1); imshow(uint8(gray_monarch)); title('Gray Monarch');
subplot(2,3,2); imshow(uint8(log_fft_mag)); title('FFT Magnitude');
subplot(2,3,3); imshow(uint8(disp_fft_phase)); title('FFT Phase');
subplot(2,3,4); imshow(uint8(ifft2D_monarch)); title('Inverse FFT');
subplot(2,3,5); imshow(disp_ifft_mag); title('IFFT of Magnitude');
subplot(2,3,6); imshow(disp_ifft_phase); title('IFFT of Phase');

% Problem 1.2
figure('Name', 'Problem 1.2 - DCT');
subplot(2,3,1); imshow(uint8(gray_monarch)); title('Gray Monarch');
subplot(2,3,2); imshow(uint8(dct_mag)); title('Normal DCT');
subplot(2,3,3); imshow(uint8(idct_monarch)); title('Normal IDCT');
subplot(2,3,4); imshow(uint8(dct16_mag)); title('DCT block 16x16');
subplot(2,3,5); imshow(uint8(idct16_monarch)); title('IDCT block 16x16');

% Problem 1.3
figure('Name', 'Problem 1.3 - DHT');
subplot(1,3,1); imshow(uint8(gray_monarch)); title('Gray Monarch');
subplot(1,3,2); imshow(uint8(dht_mag)); title('DHT Magnitude');
subplot(1,3,3); imshow(uint8(idht_monarch)); title('IDHT Monarch');

% Problem 1.5
figure('Name', 'Problem 1.5 - DWT');
subplot(1,3,1); imshow(uint8(gray_monarch)); title('Gray Monarch');
subplot(1,3,2); imshow(uint8(dwt_1)); title('DWT');
subplot(1,3,3); imshow(uint8(idwt_monarch)); title('IDWT');

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

%% Discrete Hadamard Transform (DHT)
function A = DHT(X)
  N = length(X);
  n = log2(N);

  % 2x2 Hadamard Kernel
  H1 = [1 1; 1 -1] / sqrt(2);

  % NxN Hadamard Kernel
  Hn = zeros(N,N);
  Hn(1:2,1:2) = H1;
  for i = 2:n
    j = i - 1;
    Hn(1:2^i,1:2^i) = [Hn(1:2^j,1:2^j) Hn(1:2^j,1:2^j); Hn(1:2^j,1:2^j) -Hn(1:2^j,1:2^j)] / sqrt(2);
  end
  Hn = Hn / sqrt(N);

  A = ((X*Hn).'*Hn).'*N;
end 

%% Inverse Dicrete Hadamard Transform (IDHT) Function
function A = IDHT(X)
  N = length(X);
  n = log2(N);

  % 2x2 Hadamard Kernel
  H1 = [1 1; 1 -1] / sqrt(2);

  % NxN Hadamard Kernel
  Hn = zeros(N,N);
  Hn(1:2,1:2) = H1;
  for i = 2:n
    j = i - 1;
    Hn(1:2^i,1:2^i) = [Hn(1:2^j,1:2^j) Hn(1:2^j,1:2^j); Hn(1:2^j,1:2^j) -Hn(1:2^j,1:2^j)] / sqrt(2);
  end
  Hn = Hn / sqrt(N);

  A = ((X*Hn).'*Hn).'*N;
end

%% DWT function
function A = DWT2D(X)
  [m, n] = size(X);
  half = m / 2;
  A = zeros(m,m);

  % Down sampling by 2
  for i = 1:1:half
    A1 = (X(:,2*i-1) + X(:,2*i)) / 2;
    A2 = (X(:,2*i-1) - X(:,2*i)) / 2;
    A(:,i) = A1;
    A(:,half+i) = A2;
  end
end

%% IDWT function
function A = IDWT2D(X)
  [m, n] = size(X);
  half = m / 2;
  A = zeros(m,m);

  % Up sampling by 2
  for i = 1:1:half
    A(:,2*i-1) = X(:,i) + X(:,i+half);
    A(:,2*i) = X(:,i) - X(:,i+half);
  end
end
