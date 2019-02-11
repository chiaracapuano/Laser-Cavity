% This program tries to evaluate the optical cavity transmission
% function in frequency and time space. The reference paper is:
% "Response of a ring-down cavity to an arbitrary excitation"
% Hodges, JT; Looney, P; van Zee RD.
% J. Chem. Phys., 105, 10278, (1996)
% In particular, I assume an initial gaussian excitation of the cavity modes,
% and I see what happens combining the pulse with the cavity response function (eq10)
% and (eq8). I then compare the results with the transmission function
% (eq11), which gives, anyway, less accurate results than (eq10).


clear all
format long
c=physconst('LightSpeed')*100;
nn = 2^(13); % number of steps

q_max =1;
q_in=1;
m_max=0; %considered degenerate
n_max=0;
thetax=0.02; %tilt angle for both x and y directions
thetay=0.0;
x00=0.0; %misplacement in the x and y durections
y00=0.0;

% % cavity properties
lambda=2.66*10^(-5); %cm, radiation wavelength
L_s=14; %cm, sample length
alpha_acetone=5.06*10^(-12); %1mln molecules of acetone alpha
cavity_loss = 0.002 +2*alpha_acetone*L_s ;
cavity_length = 79.3; % in cm
sample_length=14;
mirror_curvature = 100;
% step=0.5*10^(9);


RTT = 2*cavity_length/c;
f_max=(20)/RTT;
df = 2.0*f_max/(nn); % frequency spacing in Hz
t_max = 0.5/df; % maximum time given by ft limit


%Gaussian pulse characteristics
A = 1.0;
FWHM = 4*10^(-10); %Delta_time: in seconds
t0 = 0; % gaussian pulse time shift
n0 =0.3*10^9;%
w0=0.02;%cm standard deviation sigma of 0.010 cm; as w0=2*sigma=0.02 . %%used for defining which modes the cavity can sustain at 266nm
k=2*pi/lambda;
z0=pi*w0^2/lambda;

%weigths for the TEM modes that constitute the laser pulse
C(1)=1;%00
C(2)=0.2;%01
C(3)=0.2;%10
C(4)=0.1;%11



%integral of the sopace part
n=100;
xmax=0.06;
ymax=0.06;
xvar=linspace(-xmax,xmax,n);
yvar=linspace(-ymax,ymax,n);
%filestosave the distributions
fileIDx=fopen('DistributionZ00theta002.txt','w');
fileIDy=fopen('DistributionY00theta002.txt','w');

for z=-7:7
    if z>0
        thetax;
        thetay;
    else
        thetax=-thetax;
        thetay=-thetay;
    end
    
    z
    
    %function [ prodotto_space_part_output,prod_ift] = SpaceandTimePart( gaussian_pulse_ft ,dt,cavity_loss,cavity_length,mirror_curvature,lambda,w0,x00,y00,z0,z0x,z0y,thetax,thetay,z,m_max,n_max,q_in,q_max,xvar,yvar,C,var,saveFile,fileIDx,fileIDy)
    
    [out_space,out_time]=SpaceandTimePart( A,FWHM,t0,n0,cavity_loss,cavity_length,mirror_curvature,lambda,w0,x00,y00,z0,0,0,thetax,thetay,abs(z),m_max,n_max,q_in,q_max,xvar,yvar,C,t_max,nn,true,fileIDx,fileIDy);
    imagesc(xvar,yvar,out_space);
    
    xlabel('x (cm)');
    ylabel('y (cm)');
    hold on
    
    integral=sum(sum(out_space))*(ymax*2/n)^2
    
end
