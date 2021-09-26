clearvars;  close all; clc;
N = 512;
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


gaussian_beam = E0*(w0./w).*exp(-(X.^2 + Y.^2)/w^2).*exp(1i*(k*z - atan(z/zr) + k*(X.^2 + Y.^2)/(2*R)));
plane_wave = 4*(min(min(abs(gaussian_beam)))+ max(max(abs(gaussian_beam)))); %ones(N,N)*exp(i*k*pi*0.15);
 I = abs(gaussian_beam + plane_wave).^2;

figure(1);
subplot(2,2,1);
imshow(abs(gaussian_beam).^2,[]);
title('Intensity');
subplot(2,2,2);
a1 = angle(gaussian_beam);
a1 = mat2gray(a1);
imshow(a1);
title('Phase');
subplot(2,2,3);
imshow(I,[]);
title('Holographic Image');
subplot(2,2,4);
a2 = angle(gaussian_beam + plane_wave);
a2 = mat2gray(a2);
imshow(a2);
title('Diffraction Pattern');

figure(2);
subplot(1,2,1);
imshow(I,[]);
title('(a)');
subplot(1,2,2);
imshow(a2,[]);
title('(b)');



 I = I./abs(plane_wave).^2 - 1;
% %I = I - plane_wave.*conj(gaussian_beam);
 I = (I).*conj(plane_wave);

prop = Propagator(N,lambda,side_length,0.18);
prop2 = Propagator(N,lambda,side_length,0.15);
%I = abs(ifft2(ifftshift(fftshift(fft2(gaussian_beam)).*prop))).^2;

I_fft = fft2(I);
I_fft = fftshift(I_fft);

rec1 = ifft2(ifftshift(I_fft.*(prop)));
rec1 = abs(rec1).^2;
rec_fft = fftshift(fft2(rec1));
 rec = ifft2(ifftshift(rec_fft.*conj(prop2)));
 rec_int = mat2gray(abs(rec).^2); 

figure(3);
subplot(1,2,1);
imshow((rec_int),[]);
title('Intensity');
a3 = angle(rec);
a3 = mat2gray(a3);
subplot(1,2,2);
imshow(a3,[]);
title('Phase');
% 
% 
% % figure(4);
% % subplot(1,2,1);
% % imshow(abs(rec_1).^2,[]);
% % title('Intensity');
% % a4 = angle(rec_1);
% % a4 = mat2gray(a4);
% % subplot(1,2,2);
% % imshow(a4);
% % title('Phase');
% 


function [out] = FT2Dc(in)
[Nx Ny] = size(in);
f1 = zeros(Nx,Ny);
for ii = 1:Nx
    for jj = 1:Ny
         f1(ii,jj) = exp(i*pi*(ii + jj));
    end
end
FT = fft2(f1.*in);
out = f1.*FT;

end

function [out] = IFT2Dc(in)
[Nx Ny] = size(in);
f1 = zeros(Nx,Ny);
for ii = 1:Nx
    for jj = 1:Ny
        f1(ii, jj) = exp(-i*pi*(ii + jj));
    end
end
FT = ifft2(f1.*in);
out = f1.*FT;
end




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
