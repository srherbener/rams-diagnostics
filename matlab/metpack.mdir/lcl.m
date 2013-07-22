% calculate air parcel according to Clausius-Claperon Equation
%
% define initial values
p0 = 100.; %surface pressure in kPa
z0 = 0.; %height
t0 = 28+273.15; %Temperature at surface
q0 = 0.015; %1.7e-02; %absolute moisture units in kg/kg (moist air)
lclh=0; %initial value

vlev    = 31;
dp      = 1.;

%constants
g = 9.81;                %m/s2
p0 = 100.00;             %kPa
epsilon = 0.622;         %g water / g air = Rd/Rv
Tabs = 273.0;            %K
cp = 1004.67;            %J/ kg K
Rd = 287.053;            %J/ kg K
Rv = 461.50;             %J/ kg K
Lv = 2.5e6;              %J/kg
e0 = 0.611;              %kPa
gd = -(g/cp);        %K/m

for k=1:vlev-1
 p(k)     = p0 - (k-1)*dp;
 theta(k) = t0; 
 exn(k)   = (p0/p(k))^(Rd/cp);
 t(k)     = theta(k)/exn(k);
 rho(k)   = p(k)/(Rd*t(k));
 dz(k)    = dp/(g*rho(k));

 if (k==1)
 z(k) = z0;
 else 
 z(k)=z(k-1) + dz(k); 
 end

 q(k)     = q0;
 es(k)    = e0 * exp((Lv/Rv)*((1/Tabs)-(1/t(k)))); %saturation water vapor pressure
 qs(k)    = (epsilon*es(k))/p(k); %saturation mixing ratio
 gm(k)    = -gd*((1+((qs(k)*Lv)/(Rd*t(k))))/... 
	         (1+(((Lv^2)*qs(k))/(cp*Rv*(t(k)^2)))));
 end
 
%  parcel
for k=1:vlev-1
 if (z(k) <= lclh)
 thetap(k)  = theta(k);
 else
 thetap(k)  = theta(k) + gm(k)*(z(k)-lclh);
 end
end

%calculate lclh ~ height of PBL for clear skies

for k=vlev-1:-1:1
 if (q(k) >= qs(k))
 lclh = z(k);
 qhh=q(k);
 end
end
figure; plot(q,z,'r',qs,z,'--k'); hold on; plot(qhh,lclh,'o','MarkerSize', 20, 'LineWidth',2)
legend('absolute moisture at surface','saturation mixing ratio')
xlabel('moisture [kg/kg]'); ylabel('height above ground [m]')