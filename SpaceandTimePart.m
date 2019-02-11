function [ prodotto_space_part_output,prod_ift] = SpaceandTimePart( A,FWHM,t0,n0,cavity_loss,cavity_length,mirror_curvature,lambda,w0,x00,y00,z0,z0x,z0y,thetax,thetay,z,m_max,n_max,q_in,q_max,xvar,yvar,C,t_max,nn,saveFile,fileIDx,fileIDy)

c=physconst('LightSpeed')*100;

[wl,wlx,wly,etax,etay,eta, Gamma ,cosdt]=parametersSpace(cavity_loss,cavity_length,mirror_curvature,lambda,w0,x00,y00,z0,z0x,z0y,thetax,thetay,z );
k=2*pi/lambda;
[dt,t,var] = parametersTime( t_max,nn );

    
[gaussian_pulse, gaussian_pulse_ft] = GaussianPulse(A,FWHM,t0,n0,dt,t);
    
    
    
    
    
    
%x and y part of TEM modes
HG_polyx=@(m,x)   1/sqrt((2^m*factorial(m)*sqrt(pi)*wlx/sqrt(2)))*hermiteH(m,sqrt(2)*(x*cosd(thetax)-(z)*sind(thetax)- x00)/wlx).*exp(-(x*cosd(thetax)-abs(z)*sind(thetax)- x00).^2/wlx.^2).*exp(1i*(m+1/2)*etax).*exp(-1i/2*k*abs((z)*cosd(thetax)*cosd(thetay)));
HG_polyy=@(n,y)   1/sqrt((2^n*factorial(n)*sqrt(pi)*wly/sqrt(2)))*hermiteH(n,sqrt(2)*(y*cosd(thetay)-(z)*sind(thetay)- y00)/wly).*exp(-(y*cosd(thetay)-abs(z)*sind(thetay)- y00).^2/wly.^2).*exp(1i*(n+1/2)*etay).*exp(-1i/2*k*abs((z)*cosd(thetax)*cosd(thetay)));

%x and y part of laser beam
Laserbeamx=@(i,x,m) sqrt(sqrt(C(i))).*1/sqrt((2^m*factorial(m)*sqrt(pi)*wl/sqrt(2))).*hermiteH(m,sqrt(2)*(x)/wl).*exp(-(x).^2/wl.^2).*exp(1i*(m+1/2)*eta).*exp(-1i/2*k*abs((z)));
Laserbeamy=@(i,y,n) sqrt(sqrt(C(i))).*1/sqrt((2^n*factorial(n)*sqrt(pi)*wl/sqrt(2))).*hermiteH(n,sqrt(2)*(y)/wl).*exp(-(y).^2/wl.^2).*exp(1i*(n+1/2)*eta).*exp(-1i/2*k*abs((z)));

%x and y part of the space coupling coefficients
to_intx=@(m,x,i) HG_polyy(m,x).*  Laserbeamx(i,x,m);
to_inty=@(n,y,i) HG_polyy(n,y).*Laserbeamy(i,y,n);
f=@(n,i)  integral(@(y) to_inty(n,y,i), -Inf,Inf);
g=@(m,i)  integral(@(x) to_intx(m,x,i), -Inf,Inf);



prodotto=0;
prodotto_space_part=@(x,y) 0;
prodotto_space_partX=@(x) 0;
prodotto_space_partY=@(y) 0;
frequencies=[];
i=0;
idx=0;

for n=0:m_max
    boolVar=false;
    for j=0:n_max
        if(boolVar==true && n<m_max)
            break
        end
        if(boolVar==true && n<m_max)
            break
        end
        i=i+1;
        for Q = q_in:q_max
            idx=idx+1;
            frequenza=((pi*c)/cavity_length*(Q+2/pi*cosdt*(abs(j)+abs(n)+1)));
            frequencies(idx)=frequenza;
            indic=sum(frequencies(:,:)==frequenza);
            %                     if(indic>1 && j<n_max) % uncomment for not counting the degenerate modes just once (01 and 10 ie...)
            %                         boolVar=true;
            %                         break
            %                     end
            prodotto= prodotto+H(Gamma,cosdt,cavity_length,cavity_loss,var,Q,n,j);               
            prodotto_space_part=@(x,y) prodotto_space_part(x,y)+(g(n,i).*f(j,i).*HG_polyx(n,x)'.*HG_polyy(j,y));                 
            prodotto_space_partX=@(x)prodotto_space_partX(x)+ g(n,i).*HG_polyx(n,x);
            prodotto_space_partY=@(y)prodotto_space_partY(y)+f(j,i).*HG_polyy(j,y);

        end
    end
end
sommaa= gaussian_pulse_ft.*prodotto;
prod_ift = fftshift(ifft(ifftshift(sommaa)))/dt;
prodotto_space_part_output=conj(prodotto_space_part(xvar,yvar)).*prodotto_space_part(xvar,yvar);
if (saveFile==true)
       fprintf(fileIDx,'%d \n', conj(prodotto_space_partX(xvar)).*prodotto_space_partX(xvar)/sum(conj(prodotto_space_partX(xvar)).*prodotto_space_partX(xvar)));
       fprintf(fileIDy,'%d  \n',conj(prodotto_space_partY(yvar)).*prodotto_space_partY(yvar)/sum(conj(prodotto_space_partY(yvar)).*prodotto_space_partY(yvar)));

end
end

