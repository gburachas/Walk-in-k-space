%
%   EPIwalk
%
% (c) 2007-2011 Giedrius T. Burachas

clear;
%close all;
figno = 1;
%figure(figno+2); clf;

stepwise = 1;
simulatePhaseImages =1;
reconstruct =0;

% initial gradient
nlines = 8 ;% 32
starter = [zeros(2,5) [-4*ones(1,25); -4*ones(1,25)] zeros(2,5)]; 
% 2 adjacent k-space lines
blipamp = 0.3125* 32/nlines;
k2line = [2 * [ones(1,5)*blipamp zeros(1, 64) ones(1,5)*blipamp zeros(1, 64)]; [zeros(1,5) ones(1,64)*3.125 zeros(1,5) -1*ones(1,64)*3.125]]; 

% gradient timecouse

g = [starter repmat(k2line,1, nlines/2)]' * 0.75; % repmat should be 32 times (8 for demo)
b1 = zeros(length(g),1);

%subplot(3,1,3); plot(g+ones(1,length(g))'*[5 -5]); axis([0 length(g) -20 10]);
%subplot(3,1,1); plot(abs(b1));   axis([0 length(b1) min(abs(b1))*1.1-0.1 max(abs(b1))*1.1+0.1]);
%subplot(3,1,2); plot(phase(b1)); axis([0 length(b1) -pi*1.1 pi*1.1]);



%g = [zeros(size(g)); g];
%g = [g; g]';

%b1 = 1.5708*b1/max(b1)/4;

%x = [-4:.3:4];
x = [-9.5:.3:9.5];
nx = length(x);
[xx,yy] = meshgrid(x,x);
%y =x;1
y = repmat(x,length(x),1); y=y(:)';
x = [repmat(x,1,length(x)); y]';

% phantom
phRadius = 6;
phInRadius = 3;
indCirc = zeros(size(y));
indx = (x(:,1).^2+x(:,2).^2 < phRadius^2) & (x(:,1).^2+x(:,2).^2 > phInRadius^2);
indCirc(indx) = 1;


T = 0.000004; % 0.000004;% sec
T1 = 1;% sec
T2 = 0.1; % sec

%f = [-250:5:250];
f = [-100,0,100]; %[-100:10:100]; %[-10:5:10];
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
    [mx,my,mz,mrsignal] = signalbloch(b1,g,t,T1,T2,f,x,0,mx,my,mz); 
    toc
    mxy=mx+i*my;
    
    
    subplot(2,2,1);
    imagesc(reshape(real(mxy(:,ceil(length(f)/2))),XYlen,XYlen));    axis equal; axis off; colorbar; 
    subplot(2,2,2);
    imagesc(reshape(imag(mxy(:,ceil(length(f)/2))),XYlen,XYlen));    axis equal; axis off; colorbar; 
    
   
    figure(3)
    xrange =[length(starter)+1 length(mrsignal)];
    %xrange =[1 length(starter)];
    %xrange = length(starter)+[1 length(k2line)];    
    
    subplot(3,1,1);
    plot(real(mrsignal)); 
    axis([xrange min(abs(mrsignal([xrange(1):xrange(2)]))) max(real(mrsignal([xrange(1):xrange(2)])))]);
    
    subplot(3,1,2);
    plot(imag(mrsignal)); 
    axis([xrange min(phase(mrsignal([xrange(1):xrange(2)]))) max(imag(mrsignal([xrange(1):xrange(2)])))]);
    
    subplot(3,1,3);
    plot(g+ones(1,length(g))'*[5 -5]); 
    axis([xrange -20 10]);   
    
else
    
    [mx,my,mz] = signalbloch(b1(1),g(1,:),t,T1,T2,f,x,0,mx,my,mz); 

    mxy = zeros(length(x),length(f),length(g));
    mxy(:,:,1)=mx+i*my;

    mrsignal = zeros(size(b1));
    
    tic
    dur=1;
    
    for n=(dur+1):dur:(length(g)-dur+1)
        [mx,my,mz] = signalbloch(b1(n:(n+dur-1)),g(n:(n+dur-1),:),t,T1,T2,f,x,0,mx,my,mz); 
        mxy(:,:,n)=mx+sqrt(-1)*my;
        
        mrsignal(n) = sum(sum(mx)) + sqrt(-1) * sum(sum(my));
        
    end;
    toc

    if ~simulatePhaseImages
    subplot(2,2,1); % mx
    foo = reshape(squeeze(real(mxy(:,ceil(length(f)/2),n))),XYlen,XYlen);
    imagesc(foo,[-1 1]); title('X'); colorbar; axis equal; axis off; %drawnow;
    
    subplot(2,2,2); % mx
    imagesc(fftshift(abs(fft2(foo)))); title('FFT(X)'); colorbar; axis equal; axis off; %drawnow;
    
    subplot(2,2,3); % mx
    foo = reshape(squeeze(imag(mxy(:,ceil(length(f)/2),n))),XYlen,XYlen);
    imagesc(foo,[-1 1]); title('Y'); colorbar; axis equal; axis off; %drawnow;
    end;

    figure(3)
    xrange =[length(starter)+1 length(mrsignal)];
    %xrange =[1 length(starter)];
    %xrange = length(starter)+[1 length(k2line)];    
    
    subplot(2,1,1);
    plot(abs(mrsignal)); 
    axis([xrange min(abs(mrsignal([xrange(1):xrange(2)]))) max(real(mrsignal([xrange(1):xrange(2)])))]);
    
    %subplot(3,1,2);
    %plot(imag(mrsignal)); 
    %axis([xrange min(phase(mrsignal([xrange(1):xrange(2)]))) max(imag(mrsignal([xrange(1):xrange(2)])))]);
    
    subplot(2,1,2);
    plot(g+ones(1,length(g))'*[5 -5]); 
    axis([xrange -20 10]);   

% 

if reconstruct
  gmax = max(g(:,2));
  
  % sampling regions
  sampind = find(abs(g(:,2)) == gmax);
  recImage = abs(fftshift(ifft2(reshape(mrsignal(sampind),XYlen,XYlen))));
  
  figure(figno); subplot(2,2,4);
  imagesc(recImage); colorbar; axis equal; axis off; 
  title('Reconstructed image')
  
end; % reconstruct


if simulatePhaseImages % shows evolution of magnetization
    kspace = zeros(XYlen,XYlen);

for im=1:dur:length(g)
    
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
end % simulatePhaseImages

end; % stepwise

