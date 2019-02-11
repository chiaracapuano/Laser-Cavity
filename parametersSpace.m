function [ wl,wlx,wly,etax,etay,eta, Gamma ,cosdt] = parametersSpace(cavity_loss,cavity_length,mirror_curvature,lambda,w0,x00,y00,z0,z0x,z0y,thetax,thetay,z,t_max,nn )
c=physconst('LightSpeed')*100;
zx0=(+x00+(z0)*cosd(thetax)*cosd(thetay));%*sin(thetax)+z*cos(thetax)*cos(thetay);
zy0=(+y00+(z0)*cosd(thetax)*cosd(thetay));
w0x=sqrt(lambda*zx0/pi);
w0y=sqrt(lambda*zy0/pi);

wl= w0*sqrt(1+((abs(z))./(z0)).^2);%waist function at "z" (where I am calculating everything  else)
wlx=w0x*sqrt(1+(((abs(z)*cosd(thetax)*cosd(thetay))-z0x)./zx0).^2);
wly= w0y*sqrt(1+(((abs(z)*cosd(thetax)*cosd(thetay))-z0y)./zy0).^2);

etax=atan(((abs(z)*cosd(thetax)*cosd(thetay))-z0x)./zx0);
etay= atan(((abs(z)*cosd(thetax)*cosd(thetay))-z0y)./zy0);
eta= atan((abs(z))/(z0));
Gamma = -log(1-cavity_loss)/(2*cavity_length/c);
cosdt = atan(cavity_length/(sqrt(cavity_length*(2*mirror_curvature-cavity_length))));
end

