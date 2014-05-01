%
% SpiralWalk

% (c) 2007-2011 Giedrius T. Burachas

clear;
%close all;
figno = 1;
%figure(figno+2); clf;

stepwise =1;

b1 = zeros(1,1175)';

%subplot(3,1,1); plot(abs(b1));
%subplot(3,1,2); plot(phase(b1)); axis([0 length(b1) -pi*1.1 pi*1.1])
starter = []; %zeros(2,25); %[-4*ones(1,25); -4*ones(1,25)];
t=1:1165;
sprl = archiSpiral(t)/100;
g = [zeros(2,5) starter zeros(2,5) [real(sprl); imag(sprl)]]' * 0.15; % repmat should be 32 times
%subplot(3,1,3); plot(g+ones(1,length(g))'*[5 -5]);

%g = [zeros(size(g)); g];
%g = [g; g]';

%b1 = 1.5708*b1/max(b1)/4;SpiralWalk


%x = [-4:.3:4];
x = [-9.5:.3:9.5];
nx = length(x);
[xx,yy] = meshgrid(x,x);

%y =x;
y = repmat(x,length(x),1); y=y(:)';
x = [repmat(x,1,length(x)); y]';

% phantom
phRadius = 5;
indCirc = zeros(size(y));
indx = x(:,1).^2+x(:,2).^2 < phRadius^2;
indCirc(indx) = 1;

T = 0.000004; % 0.000004;% sec
T1 = 1;% sec
T2 = 0.08; % sec

%f = [-250:5:250];
f = 0; %[-100,0,100]; %[-100:10:100]; %[-10:5:10];
t = T; %[1:length(b1)]*T;

XYlen = fix(sqrt(length(x)));

mx = zeros(length(x),length(f));
my = mx;
mz = mx;
%mz = repmat(indCirc',1,length(f));
mx = repmat(indCirc',1,length(f));


figure(figno); clf; colormap(gray);

if ~stepwise
    tic    
    [mx,my,mz,mrsig] = signalbloch(b1,g,t,T1,T2,f,x,0,mx,my,mz); 
    toc
    mxy=mx+i*my;
    
    imagesc(abs(mxy));    
else

    [mx,my,mz] = signalbloch(b1(1),g(1),t,T1,T2,f,x,0,mx,my,mz);

    mxy = zeros(length(x),length(f),length(g));
    mxy(:,:,1)=mx+i*my;5

    tic
    for n=2:length(g)
        [mx,my,mz] = signalbloch(b1(n),g(n,:),t,T1,T2,f,x,0,mx,my,mz); 
        mxy(:,:,n)=mx+i*my;
    end;
    toc


%imagesc(squeeze(abs(mxy(:,ceil(length(f)/2),:))));

kspace = zeros(XYlen,XYlen);


for im=1:5:length(g)
    
    %subplot(2,2,1); % abs
    %foo = reshape(squeeze(abs(mxy(:,ceil(length(f)/2),im))),XYlen,XYlen);
    %imagesc(foo,[-1 1]); title(['ABS (t=',num2str(im),')']); colorbar; %drawnow;
    %imagesc(fftshift(abs(fft2(foo))));
    
    
    subplot(2,2,3); % mx
    fooAbs = reshape(squeeze(real(mxy(:,ceil(length(f)/2),im))),XYlen,XYlen);
    imagesc(fooAbs,[-1 1]); title('X component'); axis equal; axis off; % colorbar; %drawnow;
    
    subplot(2,2,4); % mx
    kspace = max(kspace, fftshift(abs(fft2(fooAbs))));
    imagesc(kspace, [0 1000]); title('FFT(X)'); axis equal; axis off;%colorbar; %drawnow;
    %imagesc(fftshift(abs(fft2(fooAbs)))); title('FFT(X)'); axis equal; axis off;% colorbar; %drawnow;
    title('K space sample')
    axis off; axis equal;

    subplot(2,2,2); % mx
    fooRe = reshape(squeeze(real(mxy(:,ceil(length(f)/2),im))),XYlen,XYlen);
    fooIm = reshape(squeeze(imag(mxy(:,ceil(length(f)/2),im))),XYlen,XYlen);
    %imagesc(foo,[-1 1]); title('Y'); axis equal; axis off; % colorbar; %drawnow;
    quiver(xx(nx/4:nx/2,nx/4:nx/2),yy(nx/4:nx/2,nx/4:nx/2),fooRe(nx/4:nx/2,nx/4:nx/2),fooIm(nx/4:nx/2,nx/4:nx/2)), hold off, axis image
    title('Spin phase')
    axis off; axis equal;
    
    subplot(2,2,1); % mx
    plot(g+ones(1,length(g))'*[2 -2]); hold on; axis([0 length(g) -5 5]);
    plot([im im], [-5 5]); hold off;
    %foo = reshape(squeeze(phase(mxy(:,ceil(length(f)/2),im))),XYlen,XYlen);
    %imagesc(foo,[-pi pi]); title('PHASE'); colorbar;
    title('Gradient waveforms')
    xlabel('Time (s)');ylabel('Gradient strength');
    axis off;
    drawnow;

end;


end; % stepwise

% plot of excitation at 0 freq offset


