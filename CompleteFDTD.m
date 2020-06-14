%% 2D FDTD code with unidirectional soft source and absorbing PML boundary conditions. This code can be used to simulate optical phenomena.

% Users can adjust the source wavelength, as well as the object geometry
% and material properties. 
%        
%           
%           -------------------------------------------------------
%           |                    BACK PML                         |
%           -------------------------------------------------------
%           |          |                              |           |
%           |          |                              |           |
%           | Left PML |         MAIN GRID            | Right PML |
%           |          |                              |           |
%           |          |                              |           |
%           -------------------------------------------------------
%           |                    FRONT PML                        |
%           -------------------------------------------------------   

% Property of Ana Ciocoiu, University of British Columbia
% Feb 10 2020


close all
clear all
                                                                                              
%% 1: Fundamental constants
epsilon_0 = 8.85e-12;                % permittivity of free space    
mu_0 = 4*pi*1e-7;                    % permeability of free space
c = 1/sqrt(epsilon_0*mu_0);          % speed of light (m/s)
eta_0 = sqrt(mu_0/epsilon_0);        % impedance of free space (Ohms) 
 
mode= input("Please select simulator mode: 1 is TE, 0 is TM: >>")                             % select simulator mode (1 is TE, 0 is TM)
 
%% 2: Units
micrometers = 1e-6;                  % the units used will affect the grid spacing as well as the PML boundary conditions. 
 
 
%% 3: Grid parameters
% This section defines the spacing of the 2-D grid, and sets the time step
% increments to satisfy stability (Courant) conditions.
%
delta_x = 6*micrometers;                    % spatial step increment
delta_z = delta_x;                          %
delta_y= delta_z;                           
delta_t = delta_z/(2*c);                    % Courant condition for stabiltiy of simulator
Nz = 500;                                   % spatial grid size
Nx = 500;
Nt = 3000;                                 % number of iterations in time. Nt*delta_t = total duration (seconds)
 
x = [1:Nx] * delta_x;                       % define space and time step arrays
y = [1:Nx] * delta_y;
z = [1:Nz] * delta_z;
t = [1:Nt] * delta_t;
 
if(mode==0)                                 % Set up 2-D grid for TE or TM modes
[Z X] = meshgrid(z,x);
else if(mode==1)
[Y X]= meshgrid (x,y);
    end
end
 
 
 
Ey_test= zeros(Nt,Nz);                     % initialize a set of test points
Hz_test=zeros(Nt,Nz);
E_array=zeros(Nx,Nz);                      % intialize blank array values for picking our E-field maxima
E_array_new=zeros(Nx,Nz);
%% 4: Source parameters. 
% This section sets up the properties of the electromagnetic source
 
source_z =60;                              % location of source as a function of location within the grid points (Nz)
wavelength = 300*micrometers;              % source central wavelength
omega = 2 * pi * c/wavelength;             % source central frequency
beam_waist = 500*micrometers;              % width of gaussian beam
beam_centre = Nx/2*delta_x;                % location of beam center
 
 
%% 5: Field initialization.
% This section sets up split-field electric and magnetic fields for a the PML boundary (absorbing layer). 

 
% TM with split-field for PML %
Hx(1:Nx,1:Nz) = 0.0;
Hz(1:Nx,1:Nz) = 0.0;
Eyz(1:Nx,1:Nz) = 0.0;
Eyx(1:Nx,1:Nz) = 0.0;
 
% TE with split-field for PML%
Ex(1:Nx,1:Nz) = 0.0;
Ey(1:Nx,1:Nz) = 0.0;
Hzx(1:Nx,1:Nz)=0;
Hzy(1:Nx,1:Nz)=0;
 
epsilon_r(1:Nx,1:Nz) = 1.0;   % relative  permittivity set to default of 1 across the grid (air)
 
mu_r(1:Nx,1:Nz) = 1.0;        % relative permeability set to default of 1 across the grid
 
sigma_x(1:Nx,1:Nz) = 0;       % initialize conductivity in both x and z directions (S/m) within the grid
sigma_mx(1:Nx,1:Nz) = 0;
sigma_z(1:Nx,1:Nz) = 0;
sigma_mz(1:Nx,1:Nz) = 0;
 
 
%% PML region
 
    PML_width=40;         % thickness of PML region around the grid. This can be increased to increase absorbancy of the PML
    PML_cond=80;          % condition for stable and functional PML. If the source frequency increases, this scales up by the same order of magnitude.
                          % You MUST check to see if any reflections occur off the
                          % PML region before running any new simulations
                          % if the default wavelength has been modified. 
   
  % This section sets PML conditions for total wave attenuation within the
  % region (derived in Computational Electrodynamics, Taflove) 
    sigma_x(1:PML_width,:) = PML_cond;                       
    sigma_x(Nx-PML_width:Nx,:) = PML_cond;
    sigma_mx(1:PML_width,:) = PML_cond*mu_0/epsilon_0;
    sigma_mx(Nx-PML_width:Nx,:) = PML_cond*mu_0/epsilon_0;
    
    sigma_z(:,1:PML_width) = PML_cond;
    sigma_z(:,Nz-PML_width:Nz) = PML_cond;
    sigma_mz(:,1:PML_width) = PML_cond*mu_0/epsilon_0;
    sigma_mz(:,Nz-PML_width:Nz) = PML_cond*mu_0/epsilon_0;
 
 
%% 6: Object definition.
%This section sets up the geometry of the simulation. 


% Currently set up to simulate a 300um diameter sphere of permittivity
% 2.25. 
front_sphere_radius=(25*delta_x); 

% define permittivity of object as a function of its geometry and location
% on the grid, depending on TE or TM mode
if mode ==0
epsilon_r( (X-Nx/2*delta_x).^2 + (Z-Nz/2*delta_z).^2 <=(front_sphere_radius).^2) = 2.25;

% define conductivity of object (in both x and z directions) as a function of its geometry and location
% on the grid for TE mo
% % sigma_x( (X-Nx/2*delta_x).^2 + (Z-Nz/2*delta_z).^2 <=(front_sphere_radius).^2 )=5e5;
% % sigma_z( (X-Nx/2*delta_x).^2 + (Z-Nz/2*delta_z).^2 <=(front_sphere_radius).^2 )=5e5;

% surface plot of object geometry for visualization purposes
surf(X*1e3,Z*1e3,epsilon_r), shading flat, view(2),xlabel('x [um]'), ylabel('z [um]'),axis square, drawnow
end

if mode ==1
epsilon_r( (X-Nx/2*delta_x).^2 + (Y-Nz/2*delta_z).^2 <=(front_sphere_radius).^2) = 2.25;
% define conductivity of object (in both x and z directions) as a function of its geometry and location
% on the grid for TE mo
% % sigma_x( (X-Nx/2*delta_x).^2 + (Z-Nz/2*delta_z).^2 <=(front_sphere_radius).^2 )=5e5;
% % sigma_y( (Y-Nz/2*delta_x).^2 + (Z-Nz/2*delta_z).^2 <=(front_sphere_radius).^2 )=5e5;

% surface plot of object geometry for visualization purposes
surf(X*1e3,Y*1e3,epsilon_r), shading flat, view(2),xlabel('x [um]'), ylabel('y [um]'),axis square, drawnow
end 





%% 7: FDTD factors.
%  This section is setting up mulitpliers for fields of the FDTD kernel. These 'factors' are obtained by discretizing Maxwell's equations 
% in a manner according to the FDTD algorithm. For more information, see
% Computational Electrodynamics by Taflove. 


Hxfactor1 = 1 - sigma_mx * delta_t./(2*mu_r*mu_0);
Hxfactor2 = 1 + sigma_mx * delta_t./(2*mu_r*mu_0);
Hzfactor1 = 1 - sigma_mz * delta_t./(2*mu_r*mu_0); 
Hzfactor2 = 1 + sigma_mz * delta_t./(2*mu_r*mu_0);

Exfactor1 = 1 - sigma_x * delta_t./(2*epsilon_r*epsilon_0);
Exfactor2 = 1 + sigma_x * delta_t./(2*epsilon_r*epsilon_0);
Ezfactor1 = 1 - sigma_z * delta_t./(2*epsilon_r*epsilon_0);
Ezfactor2 = 1 + sigma_z * delta_t./(2*epsilon_r*epsilon_0);

reduction = 10;   % for visualization purposes: this will plot every 10th step of the real-time simulation (line 223). 
figure

%% 8: FDTD kernel for TE of TM modes
if (mode== 0)
    for n = 1:Nt
                
        Eyz(2:Nx,2:Nz) = Eyz(2:Nx,2:Nz) .*Ezfactor1(2:Nx,2:Nz)./Ezfactor2(2:Nx,2:Nz) + ...
            delta_t./(epsilon_0 * epsilon_r(2:Nx,2:Nz).*Ezfactor2(2:Nx,2:Nz)).*...
            ( (Hx(2:Nx,2:Nz) - Hx(2:Nx,1:Nz-1) )/delta_z);
        
        Eyx(2:Nx,2:Nz) = Eyx(2:Nx,2:Nz) .*Exfactor1(2:Nx,2:Nz)./Exfactor2(2:Nx,2:Nz) + ...
            delta_t./(epsilon_0 * epsilon_r(2:Nx,2:Nz).*Exfactor2(2:Nx,2:Nz)).* (- (Hz(2:Nx,2:Nz) - Hz(1:Nx-1,2:Nz))/delta_x);
        
        Ey(2:Nx,2:Nz)=Eyz(2:Nx,2:Nz)+Eyx(2:Nx,2:Nz);
        
        Hx(1:Nx,1:Nz-1) = Hx(1:Nx,1:Nz-1).*Hzfactor1(1:Nx,1:Nz-1)./Hzfactor2(1:Nx,1:Nz-1)+...
            delta_t./(mu_0*mu_r(1:Nx,1:Nz-1).*Hzfactor2(1:Nx,1:Nz-1)).*...
            (Ey(1:Nx,2:Nz)-Ey(1:Nx,1:Nz-1))/delta_z;    
        
        Hz(1:Nx-1,1:Nz) = Hz(1:Nx-1,1:Nz).*Hxfactor1(1:Nx-1,1:Nz)./Hxfactor2(1:Nx-1,1:Nz)-...
            delta_t./(mu_0*mu_r(1:Nx-1,1:Nz).*Hxfactor2(1:Nx-1,1:Nz)).*...
            (Ey(2:Nx,1:Nz)-Ey(1:Nx-1,1:Nz))/delta_x;
        
% Initialize a gaussian-modulated sine (from settings on in section 4) source that propagates in the z direction        
%         
%               Hx(:,source_z) = -sin(omega*t(n)) * exp(- ((x-beam_centre)/beam_waist).^2)'/eta_0+Hx (:,source_z);
%               Eyz(:,source_z) = sin(omega*t(n)) * exp(- ((x-beam_centre)/beam_waist).^2)'+Eyz(:,source_z);

% Alternatively, initialize a  gaussian-modulated sine source that propagates in the x direction

               Hz(source_z,:)= sin(omega*t(n)) * exp(- ((x-beam_centre)/beam_waist).^2)/eta_0+ Hz(source_z,:);
               Eyx(source_z,:)= sin(omega*t(n)) * exp(- ((x-beam_centre)/beam_waist).^2)+ Eyx(source_z,:);
              
% Here we create test points to calculate whether energy is conserved in the grid. This is a calibration step that must be performed
% on any new simulation. Energy conservation is calculated by creating a
% square of test points for E and H fields, and then calculating the
% Poynting vector to give total power flow through each side of the square
% (W/m). More information on how this is done follows in section X. 
% The location of these test points can be modified to increase or decrease
% the size of the test square. Note that the test square should never
% contain the source. 

                  
                Ey_test(1:Nz,n)=Ey(1:Nz,1.5*source_z);
                Hx_test(1:Nz,n)=Hx(1:Nz,1.5*source_z);
                Hz_test(1:Nz,n)=Hz(1:Nz,1.5*source_z);
                               
                Ey_test2(n,1.5*source_z:0.9*Nz)=Ey(2*PML_width,1.5*source_z:0.9*Nz);
                Hx_test2(n,1.5*source_z:0.9*Nz)=Hx(2*PML_width,1.5*source_z:0.9*Nz);
                Hz_test2(n,1.5*source_z:0.9*Nz)=Hz(2*PML_width,1.5*source_z:0.9*Nz);
                
                Ey_test3(n,1.5*source_z:0.9*Nz)=Ey(0.9*Nx,1.5*source_z:0.9*Nz);
                Hx_test3(n,1.5*source_z:0.9*Nz)=Hx(0.9*Nx,1.5*source_z:0.9*Nz);
                Hz_test3(n,1.5*source_z:0.9*Nz)=Hz(0.9*Nx,1.5*source_z:0.9*Nz);
                
                Ey_test4(1:Nz,n)=Ey(1:Nz,0.9*Nz);
                Hx_test4(1:Nz,n)=Hx(1:Nz,0.9*Nz);
                Hz_test4(1:Nz,n)=Hz(1:Nz,0.9*Nz);
               
%                
%% This section of the code creates a real-time plot of the simulation (every 10th step by default) 
%         if mod(n,reduction) ==0
%             surf(X,Z,Ey), shading flat, axis equal, xlabel('x [m]'), ylabel('z [m]'), view([0 90]), drawnow
%                 
%         end
%
%%
% The code below will store current values of the electric field across the grid, the compare them with previous values, and 
% iteratively save only the largest values at each point in the grid. These
% max values are then used to create an intensity plot for the simulation,
% showing the peaks in the electric field that are caused by the object
% material & geometry. 

if(n> 0.8*Nt)   % starts taking measurements when steady state point is reached. The steady state point depends on the simulation, and must be set by the user. 
    E_array_new(1:Nz,1:Nx)=Ey;
    comp_arr_GT=E_array_new>E_array;    % will give a matrix of 1 where current array is larger
    comp_arr_LT=E_array_new<=E_array;   % will give a matrix of 1 where current array is smaller
    E_array_new =  E_array_new.*comp_arr_GT+ E_array.*comp_arr_LT;  % overwrite with larger values and keep previous larger values
    E_array=E_array_new;                % set the current values to previous values
    
end

    end
    
    
end


if (mode == 1)
    for n = 1:Nt
       
        Hzx(1:Nx-1,1:Nz-1) = Hzx(1:Nx-1,1:Nz-1) .*Hxfactor1(1:Nx-1,1:Nz-1)./Hxfactor2(1:Nx-1,1:Nz-1) - ...
            delta_t./(mu_0 * mu_r(2:Nx,2:Nz).*Hxfactor2(1:Nx-1,1:Nz-1)).*...
            (Ey(2:Nx,1:Nz-1) - Ey(1:Nx-1,1:Nz-1) )/delta_x ;
        
        Hzy(1:Nx-1,1:Nz-1) = Hzy(1:Nx-1,1:Nz-1) .*Hzfactor1(1:Nx-1,1:Nz-1)./Hzfactor2(1:Nx-1,1:Nz-1) - ...
            delta_t./(mu_0 * mu_r(2:Nx,2:Nz).*Hzfactor2(1:Nx-1,1:Nz-1)).*(- (Ex(1:Nx-1,2:Nz) - Ex(1:Nx-1,1:Nz-1))/delta_y);
        
        Hz(1:Nx-1,1:Nz-1)= Hzy(1:Nx-1,1:Nz-1)+ Hzx(1:Nx-1,1:Nz-1);
        
        Ey(2:Nx,2:Nz) = Ey(2:Nx,2:Nz).*Exfactor1(2:Nx,2:Nz)./Exfactor2(2:Nx,2:Nz)-...
            delta_t./(epsilon_0*epsilon_r(2:Nx,2:Nz).*Exfactor2(2:Nx,2:Nz)).*...
            (Hz(2:Nx,2:Nz)-Hz(1:Nx-1,2:Nz))/delta_x;
        
        
        Ex(2:Nx,2:Nz) = Ex(2:Nx,2:Nz).*Ezfactor1(2:Nx,2:Nz)./Ezfactor2(2:Nx,2:Nz)+...
            delta_t./(epsilon_0*epsilon_r(2:Nx,2:Nz).*Ezfactor2(2:Nx,2:Nz)).*...
            (Hz(2:Nx,2:Nz)-Hz(2:Nx,1:Nz-1))/delta_y;
        
% Initialize a gaussian-modulated sine (from settings on in section 4) source that propagates in the y direction       

%         Hzy(:,source_z) = sin(omega*t(n)) * exp(- ((x-beam_centre)/beam_waist).^2)'/eta_0+Hzy (:,source_z);
%         Ex(:,source_z) = -sin(omega*t(n)) * exp(- ((x-beam_centre)/beam_waist).^2)'+Ex(:,source_z)    


% Alternatively, initialize a  gaussian-modulated sine source that propagates in the x direction
         Hzx(source_z,:)= sin(omega*t(n)) * exp(- ((x-beam_centre)/beam_waist).^2)/eta_0+ Hzx(source_z,:);
        Ey(source_z,:)= sin(omega*t(n)) * exp(- ((x-beam_centre)/beam_waist).^2)+ Ey(source_z,:);

% Here we create test points to calculate whether energy is conserved in the grid. This is a calibration step that must be performed
% on any new simulation. Energy conservation is calculated by creating a
% square of test points for E and H fields, and then calculating the
% Poynting vector to give total power flow through each side of the square
% (W/m^2). More information on how this is done follows in section X. 
% The location of these test points can be modified to increase or decrease
% the size of the test square. Note that the test square should never
% contain the source. 

                  
                Ey_test(1:Nz,n)=Ey(1:Nz,1.5*source_z);
                Ex_test(1:Nz,n)=Ex(1:Nz,1.5*source_z);
                Hz_test(1:Nz,n)=Hz(1:Nz,1.5*source_z);
                               
                Ey_test2(n,1.5*source_z:0.9*Nz)=Ey(2*PML_width,1.5*source_z:0.9*Nz);
                Ex_test2(n,1.5*source_z:0.9*Nz)=Ex(2*PML_width,1.5*source_z:0.9*Nz);
                Hz_test2(n,1.5*source_z:0.9*Nz)=Hz(2*PML_width,1.5*source_z:0.9*Nz);
                
                Ey_test3(n,1.5*source_z:0.9*Nz)=Ey(0.9*Nx,1.5*source_z:0.9*Nz);
                Ex_test3(n,1.5*source_z:0.9*Nz)=Ex(0.9*Nx,1.5*source_z:0.9*Nz);
                Hz_test3(n,1.5*source_z:0.9*Nz)=Hz(0.9*Nx,1.5*source_z:0.9*Nz);
                
                Ey_test4(1:Nz,n)=Ey(1:Nz,0.9*Nz);
                Ex_test4(1:Nz,n)=Ex(1:Nz,0.9*Nz);
                Hz_test4(1:Nz,n)=Hz(1:Nz,0.9*Nz);


        
 %% This section of the code creates a real-time plot of the simulation (every 10th step by default) 
%        if mod(n,reduction) ==0
%             surf(X,Y,Hz), shading flat, axis equal, xlabel('x'), ylabel('y'), view([0 90]), drawnow
%         end
%         
        
        %%
% The code below will store current values of the electric field across the grid, the compare them with previous values, and 
% iteratively save only the largest values at each point in the grid. These
% max values are then used to create an intensity plot for the simulation,
% showing the peaks in the electric field that are caused by the object
% material & geometry. 

  % starts taking measurements when steady state point is reached. The steady state point depends on the simulation, and must be set by the user. 
if(n> 0.8*Nt)    
  E_array_new(1:Nz,1:Nx)=Hz;
    comp_arr_GT=E_array_new>E_array;    % will give a matrix of 1 where current array is larger
    comp_arr_LT=E_array_new<=E_array;   % will give a matrix of 1 where current array is smaller
    E_array_new =  E_array_new.*comp_arr_GT+ E_array.*comp_arr_LT;  % overwrite with larger values and keep previous larger values
    E_array=E_array_new;                % set the current values to previous values
end   


    end
    
end


%% 9: Poynting vector calculations: use with a source and no objects to test whether the grid is stable. 
% This code takes the electic field values from the simulation, which were
% sampled from around a closed square box within the grid. The time-averaged
% Poynting vector (W/m^2) is calculated for each size of the box. It is then "integrated" (dotted with each box side's normal
% vector and summed over the length of the side) to give the total amount
% of power flowing of in/out of the box. Ideally, due to conservation of
% energy, this value should be smaller than the user-set threshold value to
% ensure a stable simulation. If the threshold is not met, consider
% reducing the grid spacing.

threshold=1e-5;  %% error threshold for conservation of energy calculations

if (mode==1)
   H_poynt= [zeros(1,Nz)' zeros(1,Nz)' Hz_test(1:Nz,n)];                         
   E_poynt= [ Ex_test(1:Nz,n) Ey_test(1:Nz,n) zeros(1,Nz)' ];
   S= 0.5*real(cross(hilbert(E_poynt), conj(hilbert(H_poynt))));
   
   H_poynt2= [zeros(1,201)' zeros(1,201)' Hz_test2(n,1:201)'];            
   E_poynt2= [ Ex_test2(n,1:201)' Ey_test2(n,1:201)' zeros(1,201)' ];
    S2= 0.5*real(cross(hilbert(E_poynt2), conj(hilbert(H_poynt2))));
      
   H_poynt3= [zeros(1,201)' zeros(1,201)' Hz_test3(n,1:201)'];            
   E_poynt3= [ Ex_test3(n,1:201)' Ey_test3(n,1:201)' zeros(1,201)' ];
    S3= 0.5*real(cross(hilbert(E_poynt3), conj(hilbert(H_poynt3))));
   
   H_poynt4= [zeros(1,Nz)' zeros(1,Nz)' Hz_test4(1:Nz,n)];
   E_poynt4= [ Ex_test4(1:Nz,n) Ey_test4(1:Nz,n) zeros(1,Nz)' ];
   S4= 0.5*real(cross(hilbert(E_poynt4), conj(hilbert(H_poynt4))));
end

if (mode==0)
       
   H_poynt= [Hx_test(1:Nz,n) zeros(1,Nz)' Hz_test(1:Nz,n)];
   E_poynt= [zeros(1,Nz)' Ey_test(1:Nz,n)  zeros(1,Nz)' ];
   S= 0.5*real(cross(hilbert(E_poynt), conj(hilbert(H_poynt))));
   
   H_poynt2= [Hx_test2(n,1:201)' zeros(1,201)' Hz_test2(n,1:201)'];
   E_poynt2= [zeros(1,201)' Ey_test2(n,1:201)'  zeros(1,201)' ];
   S2= 0.5*real(cross(hilbert(E_poynt2), conj(hilbert(H_poynt2))));

   H_poynt3= [Hx_test3(n,1:201)' zeros(1,201)' Hz_test3(n,1:201)'];
   E_poynt3= [zeros(1,201)' Ey_test3(n,1:201)'  zeros(1,201)' ];
   S3= 0.5*real(cross(hilbert(E_poynt3), conj(hilbert(H_poynt3))));

   
   H_poynt4= [Hx_test4(1:Nz,n) zeros(1,Nz)' Hz_test4(1:Nz,n)];
   E_poynt4= [zeros(1,Nz)' Ey_test4(1:Nz,n)  zeros(1,Nz)' ];
   S4= 0.5*real(cross(hilbert(E_poynt4), conj(hilbert(H_poynt4))));
   
   
end

sumtotal= -sum(S(3))+sum(S4(3))+sum(S3(1))-sum(S2(1));  %check conservation of energy
if(sumtotal<threshold)
    disp('stable: hooray!')
else
    disp('not stable: try reducing grid spacing! ')
end

%% 10 : This section plots the maximum E or H-field array obtained in 8 as a heat map
if (mode ==0)
arraycolor=surf(Z/delta_x,X/delta_x, E_array_new);
shading interp;
view(2);
set(arraycolor,'linestyle','none');
end

if (mode ==1)
arraycolor=surf(Y/delta_x,X/delta_x, E_array_new);   
shading interp;
view(2);
set(arraycolor,'linestyle','none');
end

