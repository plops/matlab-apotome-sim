%% 
% all lengths are in micrometer
NA=1.38; % numerical aperture
ri=1.515;% refractive index
lambda=.5; % emission vacuum wavelength
rbfp=NA/lambda; % radius of bfp for oil immersion objective
s=[256 256]; % image size, filling the full diameter of the bfp
apertureX=2*rbfp; % size of the bfp image in um
apertureY=2*rbfp;
cosTheta=newim(s)+1;                  
kz_of_kxky=cosTheta;

klen=ri/NA; % radius of ewald sphere in units of the aperture

kx=xx(s,'freq')*(apertureX);                   % normalized to ellipse
ky=yy(s,'freq')*(apertureY);                   % normalized to ellipse               
rho2=kx^2+ky^2;

tmp2=klen^2-rho2;   % k_0 here is the ratio between k_0 in pixels and aperture in pixels
aperture=tmp2>0; % only a disk of radius rbfp is contributing
cosTheta(aperture)=sqrt(tmp2(aperture));        
                    % Theta, being the angle to the optic axis
kz_of_kxky=cosTheta;% still in units of kx and ky
%% prepare phase to propagate to different slices of the stack
nz=3; % number of slices
dz=.3; % step size in um
z=xx(nz)*dz; % defocus in um
sxy=.04095; % um per pixel
vol=[s nz];
fpropmat=newim(vol,'dcomplex');
for i=0:length(z)-1
    fpropmat(:,:,i)=exp(kz_of_kxky(:,:)*(1i*klen/sxy*z(i)));
end

%%
% the orders of the following grating would hit exactly the periphery of
% the aperture in the bfp:
% grat=sin(2*pi*xx(s,'freq')*s(1)/apertureX*klen) * (rr(s,'freq')<.1);
% so put the orders halfway between center and periphery:
grat=sin(2*pi*xx(s,'freq')*s(1)*.5/apertureX*klen) * (rr(s,'freq')<.1);

% TODO: simulate grating with 3 orders, first order should hit 2/3 of bfp
% grat=abs(sin(2*pi*xx(s,'freq')*s(1)*(2/3)/apertureX*klen)) * (rr(s,'freq')<.1);

%% get intensity of coherent stack in sample space with parallel
%% illumination
coh=newim(vol,'dfloat');
kgrat=ft(grat);
for i=0:length(z)-1
    coh(:,:,i)=abs(ift(kgrat.*squeeze(fpropmat(:,:,i)).*aperture)).^2;
end
clear kgrat

%% make aperture for illumination angles half as wide as bfp
ill_aperture=(rr(s,'freq')*apertureX/klen)<.5;
% I don't want to add several thousand illumination angles: decimate
ill_factor=8;
ill_aperture_small=(rr(s/ill_factor,'freq')*apertureX/klen)<.5;

%% for testing: this statement moves the center of the diffraction pattern
%% from 128 128 to 100 100:
%ft(grat.*exp(1i*(double(kx(100,100))*xx(s)+double(ky(100,100))*yy(s))))

%% accumulate an incoherent image for ideal spatial incoherence of light
%% source
incoh=newim(vol,'dfloat');
incoh_bfp=newim(s,'dfloat');
[row col]=size(ill_aperture_small);
all=sum(ill_aperture_small);
count=0;
for j=0:row-1
    for i=0:col-1
        if ill_aperture_small(i,j)
            % add this illumination angle to the incoherent image
            shifter=exp(1i*(kx(i*ill_factor,j*ill_factor)*xx(s)+...
                ky(i*ill_factor,j*ill_factor)*yy(s)));
            kgrat=ft(grat.*shifter);
            for k=0:length(z)-1
                incoh(:,:,k)=squeeze(incoh(:,:,k))+abs(ift(...
                    kgrat(:,:).*squeeze(fpropmat(:,:,k)).*aperture(:,:))).^2;
            end
            % accumulate all patterns in the bfp as well, that is what you
            % would see if you put a paper into the bfp. if the simulation
            % doesn't use enough angles, you will see it here
            incoh_bfp(:,:)=squeeze(incoh_bfp(:,:))+abs(kgrat(:,:)).^2;
            count=count+1;
            [i j count all] % print position and how far we are into calculation
        end
    end
end
clear kgrat shifter

%% check the bfp
overlay(incoh_bfp/max(incoh_bfp)*255,~aperture)

%% for simulating what happens to the contrast with different diameters of
%% the illumination aperture I want to reuse calculations as often as
%% possible
% I want to simulate the following sizes of the illumination aperture
% relative to the bfp
% .1 .25 .4 .5 .6 .7 .9 1. 1.1
ill_sizes=[.1 .25 .4 .5 .6 .7 .9 1. 1.1 1.2 1.3 1.4 1.5 1.6 1.7];
ill_ap=newim([s/ill_factor length(ill_sizes)],'bin');
for k=0:length(ill_sizes)-1
    ill_ap(:,:,k)=rr(s/ill_factor,'freq')*apertureX/klen<ill_sizes(1+k);
end

%% accumulate several incoherent images with varying illumination aperture
incohs=newim([vol length(ill_sizes)],'dfloat');
current=newim(vol,'dfloat');
incoh_bfps=newim([s length(ill_sizes)],'dfloat');
[row col hei]=size(ill_ap);
all=sum(ill_ap(:,:,hei-1));
count=0;
for j=0:row-1
    for i=0:col-1
        if ill_ap(i,j,hei-1) % does the biggest circle contain that angle?
            % coherent image for this angle
            shifter=exp(1i*(kx(i*ill_factor,j*ill_factor)*xx(s)+...
                ky(i*ill_factor,j*ill_factor)*yy(s)));
            kgrat=ft(grat.*shifter);
            for k=0:length(z)-1
                current(:,:,k)=abs(ift(kgrat(:,:).*squeeze(fpropmat(:,:,k)).*aperture(:,:))).^2;
            end
            % add the image to all the incoherent images that need it
            for q=0:hei-1
                if ill_ap(i,j,q) % FIXME the squeeze doesn't work with nz=1 :-(
                  incohs(:,:,:,q)=squeeze(incohs(:,:,:,q))+current(:,:,:);
                  incoh_bfps(:,:,q)=squeeze(incoh_bfps(:,:,q))+abs(kgrat(:,:)).^2;
                end
            end
            count=count+1;
            [i j count all]
        end
    end
end
clear kgrat shifter current

%% normalize each image according to the number of angles it received
norm_incohs=incohs;
for q=0:hei-1
    norm_incohs(:,:,:,q)=incohs(:,:,:,q)./sum(ill_ap(:,:,q));
end

%%
contrast=ill_sizes;
for q=0:hei-1
    selection=norm_incohs(114:143,floor(s(1)/2),floor(nz/2),q);
    contrast(q+1)=max(selection)-min(selection);
end
plot(ill_sizes,contrast)
xlabel('illumination aperture / bfp aperture')
ylabel('contrast max-min of in-focus grating image')
% print(gcf,'-depsc','/dev/shm/o.eps')
print(gcf,'-dpng','~/1029/contrast.png')

% for the email I want one big image that contains 3 slices of the stack
% next to the image in the bfp for each of the illumination situations
%% overlay an bfp aperture ontop of the incoherent images from the bfp
%% plane
bfp_im=.5*incoh_bfps/max(incoh_bfps)+.5*repmat(aperture,[1 1 length(ill_sizes)]);
%% add the bfp image to the stack
combin=newim([s nz+1 length(ill_sizes)],'dfloat');
for q=0:length(ill_sizes)-1
    combin(:,:,0:nz-1,q)=squeeze(norm_incohs(:,:,:,q));
    combin(:,:,nz,q)=bfp_im(:,:,q);
end
%% shuffle around the dimensions and convert into 2d mosaic
bigcomb=reshape(permute(reshape(combin,[256 256*4 15]),[2 1 3]),[256*4 256*15]);
% save to file
writeim(bigcomb/max(bigcomb)*255,'/home/martin/1029/apotome-simul.jpg','JPEG');