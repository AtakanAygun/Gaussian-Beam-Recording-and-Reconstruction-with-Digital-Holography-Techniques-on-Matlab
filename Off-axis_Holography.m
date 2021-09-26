clearvars;  close all; clc;
N = 256;
E0 = 1;
lambda = 632.8e-9;
z = 0.15;
k = 2*pi/lambda;
w0 = 50e-6;
zr = pi*w0^2/lambda;
w = w0*sqrt(1+ (z/zr)^2);
R = z*(1 + (zr/z))^2;
side_length = 0.003;
delta = side_length/N;
x = -side_length/2 + delta/2:delta:side_length/2 - delta/2;
y = -side_length/2 + delta/2:delta:side_length/2 - delta/2;
[X, Y] = meshgrid(x,y);
alpha = 50;
k2 = 2*pi/(1e-6);
ky = k*sin(alpha);
kz = k*cos(alpha);

plane_wave = 1*10^-2*exp(1i*kz*z)*exp(1i*ky*Y);

gaussian_beam = E0*(w0./w).*exp(-(X.^2 + Y.^2)/w^2).*exp(1i*(k*z - atan(z/zr) + k*(X.^2 + Y.^2)/(2*R)));

% I = abs(gaussian_beam + plane_wave).^2;
% figure(1);
% imshow(I,[]);
% a1 = angle(gaussian_beam + plane_wave);
% figure(2);
% a1 = mat2gray(a1);
% imshow(a1);

% figure(1);
% subplot(1,2,1);
% imshow(I,[]);
% title('(a)');
% subplot(1,2,2);
% a1 = angle(gaussian_beam);
% a1 = mat2gray(a1);
% a2 = angle(gaussian_beam + plane_wave);
% a2 = mat2gray(a2);
% imshow(a2);
% title('(b)');

r = 1:5*N;
c = 1:5*N;
[C, R] = meshgrid(c, r);
g = zeros(5*N);
g(256*2+1:256*3,256*2+1:256*3) = gaussian_beam;

prop = Propagator(N*5,lambda,side_length,0.01);
prop2 = Propagator(4*N,lambda,side_length,0.02);


g_p = fftshift(fft2(g));
g_p = ifft2(ifftshift(g_p.*(prop)));
g_p = g_p(256*2+1:256*3,256*2+1:256*3);

Az = fftshift(fft2(g_p));
Az2 = zeros(4*N);
Az2(256*3/2+1:256*3/2+256,256*3/2+1:256*3/2+256) = Az;
EOf = ifft2(ifftshift(Az2));

AV = (min(min(abs(EOf))) + max(max(abs(EOf))))/2;

r2 = 1:4*N;
c2 = 1:4*N;
[C2, R2] = meshgrid(c2, r2);
Ref =  AV*exp(1i*ky*R2);%exp(1i*2*pi*(alpha)*delta/4.*(R2-2*N-1)/lambda + 1i*2*pi*(alpha)*delta/4.*(C2-2*N-1)/lambda);

IH = (EOf + Ref).*conj(EOf + Ref);
IH = IH - abs(Ref).^2;
IH = IH(257:256*3,257:256*3);

  

r3 = 1:2*N;
c3 = 1:2*N;
[C3, R3] = meshgrid(c3, r3);
THOR = ((R3-N-1).^2 + (C3-N-1).^2).^0.5;
A= THOR.*delta/4;
QP = exp(1i*pi/lambda/0.02.*(A.^2));
FTS = fftshift(fft2(fftshift(IH.*QP)));
I2 = FTS.*conj(FTS);

g_p2 = fftshift(fft2(EOf(257:256*3,257:256*3)));
  g_p2 = ifft2(ifftshift(g_p2.*QP));
   %g_p2 = g_p2(257:256*3,257:256*3);



figure(1);
subplot(1,3,1);
imshow(g_p2.*conj(g_p2),[]);

subplot(1,3,2);
imshow(IH,[]);
subplot(1,3,3);
imshow(angle(g_p2),[]);


figure(2);
subplot(1,3,1)
subplot(1,3,2);
imshow(I2,[]);
subplot(1,3,3);
imshow(angle(FTS),[]);

figure(3);
imshow(IH,[]);
title('Off-axis Hologram');

figure(4);
subplot(2,2,1);
imshow(g_p2.*conj(g_p2),[]);
title('Object Intensity Distribution');
subplot(2,2,2);
imshow(angle(g_p2),[]);
title('Object Phase Distribution');
subplot(2,2,3);
imshow(I2,[]);
title('Reconstructed Wavefield Intensity Distribution');
subplot(2,2,4);
imshow(angle(FTS),[]);
title('Reconstructed Wavefield Phase Distribution');
figure(5);
subplot(1,2,1);
imshow(g_p2.*conj(g_p2),[]);
title('Object Intensity Distribution');
subplot(1,2,2);
imshow(angle(g_p2),[]);
title('Object Phase Distribution');
figure(6);
subplot(1,2,1);
imshow(I2,[]);
title('Reconstructed Wavefield Intensity Distribution');
subplot(1,2,2);
imshow(angle(FTS),[]);
title('Reconstructed Wavefield Phase Distribution');










 




    
%     
%     H_fft = fftshift(fft2(H));
%     
%     rec = ifft2(ifftshift(H_fft.*conj(prop2)));
%     
%     figure(1);
%     subplot(1,3,1);
%     imshow(H,[]);
%     subplot(1,3,2);
%     imshow(abs(rec),[]);
%     subplot(1,3,3);
%     imshow(angle(rec),[]);
    
% I = I - abs(plane_wave).^2;
% %I = I - plane_wave.*conj(gaussian_beam);
% I = I.*plane_wave;
% 
% z2 = 0.5;
% U = X;
% V = Y ;% - z2*tan(alpha);
% 
% u = x;
% v = y + z2*tan(alpha);
% B = zeros(512,512);
% 
% B = exp(1i*k*z2)/(1i*k*z2)*exp(1i*pi*(U.^2 + V.^2)./(lambda*z2)).*fftshift(fft2(I.*exp(1i*pi*(X.^2 + Y.^2)/(lambda*z2))));
% B = (B);
% R0 = 50;
% B_fft = fftshift(fft2(B));
% B_1fft = zeros(N,N);
% 
% for ii=1:N
%     for jj=1:N
%      
%     x = ii - N/2;
%     y = jj - N/2;
%     
%     if (sqrt(x^2 + y^2) > R0) 
%         B_1fft(ii, jj) = B_fft(ii,jj); 
%     end
%     end
% end
% 
% 
% B_1 = ifft2(ifftshift(B_1fft));
% 
% figure(2);
% subplot(1,2,1);
% imshow(abs(B).^2,[]);
% subplot(1,2,2);
% imshow(mat2gray(angle(B)),[]);
% 
% figure(3);
% subplot(1,2,1);
% imshow(abs(B_1).^2,[]);
% subplot(1,2,2);
% imshow(mat2gray(angle(B_1)));
% 


% I_fft = fft2(I);
% I_fft = fftshift(I_fft);
% prop = Propagator(N,lambda,side_length,0.03);
% prop2 = Propagator(N,lambda,side_length,0.03);
% 
% rec1 = ifft2(ifftshift(I_fft.*(prop2)));
% rec1 = abs(rec1).^2;  
% rec_fft = fftshift(fft2(rec1));
%  rec = ifft2(ifftshift(rec_fft.*conj(prop2)));
% rec_int = mat2gray(abs(rec).^2);
% % 
% figure(4);
% subplot(1,2,1);
% imshow(rec_int);
% title('Intensity');
% a2 = angle(rec);
% a2 = mat2gray(a2);
% subplot(1,2,2);
% imshow(a2);
% title('Phase');





% 
function [p] =  Propagator(N,lambda,area,z)
    p = zeros(N,N);
    
    for ii = 1:N
        for jj = 1:N
            alpha = lambda*(ii -N/2 -1)/area;
            beta = lambda*(jj - N/2 -1)/area;
            if ((alpha^2 + beta^2) <= 1)
                p(ii,jj) = exp(-2*pi*1i*z*sqrt(1 - alpha^2 - beta^2)/lambda);
            end
        end
    end
end

