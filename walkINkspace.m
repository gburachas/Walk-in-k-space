%
% walkINkspace
%
% ==== Show simulation of M(time,position,freq)

% (c) 2007-2011 Giedrius T. Burachas

clear;
%close all;
figno = 1;
stepwise = 1;
b1 = 0.7*[zeros(1,150) 1*sinc([-12.4:0.1:12.5]/2.08) zeros(1,250) 1*sqrt(-1)*sinc([-12.4:0.1:12.5]/2.08) zeros(1,175)]';
b1 = zeros(size(b1));
%g = 4*[-ones(1,125) ones(1,250) -ones(1,125) -ones(1,125)  ones(1,250) -ones(1,125) zeros(1,200)];
g = [zeros(1,25) -ones(1,125) ones(1,250) -ones(1,125)  zeros(1,550)];
%g = 0.2 * [g; [zeros(1,500) -ones(1,125)  ones(1,250) -ones(1,125) zeros(1,200)]];
g = 0.5 * [[zeros(1,550) -ones(1,125)  ones(1,250) -ones(1,125) zeros(1,25)];g]';

if 0
figure(figno+2); clf;
subplot(3,1,1); plot(abs(b1));
subplot(3,1,2); plot(phase(b1)); axis([0 length(b1) -pi*1.1 pi*1.1])
subplot(3,1,3); plot(g+ones(1,length(g))'*[2 -2]);
end
%g = [zeros(size(g)); g];
%g = [g; g]';

%b1 = 1.5708*b1/max(b1)/4;


%x = [-4:.3:4];
x = [-9.5:.3:9.5]; % takes very long to simulate
nx = length(x);
[xx,yy] = meshgrid(x,x);
%y =x;
y = repmat(x,length(x),1); y=y(:)';x = [repmat(x,1,length(x)); y]';

% defining a circular phantom
phRadius = 5;
indCirc = zeros(size(y));
indx = (x(:,1).^2+x(:,2).^2 < phRadius^2); 
indCirc(indx) = 1;


T = 0.000004; % 0.000004;% sec - the timebase for signal integration
T1 = 1;% sec
T2 = 0.08; % sec

%f = [-250:5:250];
f = 0; % [-100,0,100]; %[-100:10:100]; %[-10:5:10]; % the range of off-resonance offsets
t = T; %[1:length(b1)]*T;

XYlen = fix(sqrt(length(x)));

% magnetization vector components 
mx = zeros(length(x),length(f));  % X component
my = mx; 
mz = mx;
%mz = repmat(indCirc',1,length(f));
mx = repmat(indCirc',1,length(f));


figure(figno); clf; colormap(gray);

if ~stepwise %  one-step calculation of the full evolution of magnetization
    tic    
    [mx,my,mz,mrsignal] = signalbloch(b1,g,t,T1,T2,f,x,0,mx,my,mz); 
    %[mx,my,mz] = bloch(b1,g,t,T1,T2,f,x,0,mx,my,mz); 
    toc
    mxy=mx+i*my;
    subplot(2,2,1);
    imagesc(reshape(real(mxy(:,ceil(length(f)/2))),XYlen,XYlen));    axis equal; %colorbar; 
    subplot(2,2,2);
    imagesc(reshape(imag(mxy(:,ceil(length(f)/2))),XYlen,XYlen));    axis equal; %colorbar; 
else  % stepwise evolution of magnetization
    
    [mx,my,mz] = signalbloch(b1(1),g(1,:),t,T1,T2,f,x,0,mx,my,mz); 

    mxy = zeros(length(x),length(f),length(g));
    mxy(:,:,1)=mx+i*my;

    tic
    dur=1;
    
    for n=(dur+1):dur:(length(g)-dur+1)
        [mx,my,mz] = signalbloch(b1(n:(n+dur-1)),g(n:(n+dur-1),:),t,T1,T2,f,x,0,mx,my,mz); 
        mxy(:,:,n)=mx+i*my;
    end;
    toc


% PLOTTING
for im=1:dur:length(g)
    
    %subplot(2,2,1); % abs
    %foo = reshape(squeeze(abs(mxy(:,ceil(length(f)/2),im))),XYlen,XYlen);
    %imagesc(foo,[-1 1]); title(['ABS (t=',num2str(im),')']); colorbar; %drawnow;
    %imagesc(fftshift(abs(fft2(foo))));
    
    
    subplot(2,2,3); % mx
    fooAbs = reshape(squeeze(real(mxy(:,ceil(length(f)/2),im))),XYlen,XYlen);
    imagesc(fooAbs,[-1 1]); title('X component'); axis equal; axis off; % colorbar; %drawnow;
    
    subplot(2,2,4); % mx
    imagesc(fftshift(abs(fft2(fooAbs)))); title('FFT(X)'); axis equal; axis off;% colorbar; %drawnow;
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

    if im==1, pause; end
end;


end; % stepwise

%ylabel('Magnetization');
%xlabel('Location x');

