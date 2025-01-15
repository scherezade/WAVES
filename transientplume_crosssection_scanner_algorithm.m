close all
clear all

%% Readme

%Scherezade Barquero,
% Contact info: scherezadebarquero@gmail.com, %scherezade.bbarquero@gmail.com (mbalsera@ing.uc3m.es: unavailable)
% Last update: January, 15th - 2025

   
% Paper title: Reconstruction of the transient plume cross-section of a Pulsed Plasma
% Thruster. Authors: S. Barquero, M. Merino, J. Navarro-Cavallé.
% Plasmas and Space Propulsion Research Group - EP2, Aerospace Engineering Department, University Carlos III de Madrid (Leganés, Spain)
% Plasma Sources Science and Technology  

% Paper abstract: The exhaust of a small ablative pulsed plasma thruster (PPT), fed with polytetrafluoroethylene and operated at 1000 V
% of discharge voltage and 6 μF of capacitance, is characterized by means of a novel diagnostic system to time-reconstruct the
% cross-sectional expansion of unsteady electrical plasma thrusters. This experimental technique consists of an array of electro-
% static wire probes biased at the ion saturation regime. The two-dimensional and time-dependent ion current distribution is
% reconstructed from the probe data using a variable separation algorithm. The PPT plume contains at least three distinct ion
% groups with different mean velocities, with the second one carrying the major part of the ion current. Spatially, the plume
% exhibits a single-peaked profile in the direction perpendicular to the PPT electrodes, while in the direction parallel to them
% it features two peaks and a greater divergence angle. Small spatial asymmetries are present. In particular, there is a alight
% deviation of the current towards the cathode and one of the sides of the channel.

% This script reconstructs the plasma exhaust by computing, according to the algorithm presented in the paper*, 
% 1) the ion current density distribution maps j(x,y,t). Additionally, from j(x,y,t), are also calculated:
% 2) its row- and column-wise profiles as well as 
% 3) the time-integrated charge-surface density distribution map rho(x,y) together with its corresponding cross profiles, 
% 4) and the time-series of the total current collected by the grid
% scanner.
% *) In addition, the linearized version of the algorithm, only valid for uniform cell size, is also included in the present script. 


% Paper data available at the following URL/DOI: 10.5281/zenodo.13820945:
% They include the measured data and the exact solution of the linearized version of the inverse problem.

%%
%----------------------------------------------------------------------------
%-----------------------To Modify by the U S E R ---------------------------- 
%----------------------------------------------------------------------------
%0.1: Folder paths
%----------------------------------------------------------------------------
motherfolder = "D:\"; 
ZENODOdata = "FinalZENODO_PaperGrid.mat"; % measured data downloaded from ZENODO
computedsolution_matfile = "FinalZENODO_PaperGrid_computedsolution.mat"; % to store the ion current density maps j(x,y,t)  and related calculations

%----------------------------------------------------------------------------
%0.2 FLAGS
%----------------------------------------------------------------------------

computedsolution.flags.flag1   = 0; %To compute j
computedsolution.flags.flag1a  = 1; % Nonlinear system 
computedsolution.flags.flag1a1 = 1; % nonuniform grid (similar result than using a uniform grid)
% By default, if flag1a == 0: Linearized algorithm. This exact solution is only valid for uniform grids (flag1a1 is ignored)

computedsolution.flags.flag6s= 1; %To save j, only works if flag1 is active

computedsolution.flags.flag2 = 0; % (j+ things)     To calculate from the computed j: 1) jcenterofmass and j row- and column-wise profiles at the 3 times of interest, and 2) I_tot, and Qtot from I_tot
computedsolution.flags.flag7s= 1; %To save j+ things. Only works if flag2 is active

computedsolution.flags.flag3 = 0; % (rho+ things)   To calculate from the computed j: rho,rhocenterofmass,rho row- and column-wise profiles and Qtot from rho 
computedsolution.flags.flag8s= 1; %To save rho+ things. Only works if flag3 is active

computedsolution.flags.flag4 = 0; % (theta+ things) To compute the half-angle-at-half-width instantaneous and average divergence angles and the Qstripe/Qtot x,y profiles
computedsolution.flags.flag9s= 1; %To save theta+ things. Only works if flag4 is active

computedsolution.flags.flag10 = 1; %Postprocess on
computedsolution.flags.flag00 = 1; %For checking: to reconstruct the measured data from computed j 
computedsolution.flags.flag11s= 1; %To save the figures

%----------------------------------------------------------------------------
% 0.3: Times of interest
%----------------------------------------------------------------------------
%ION COLLECTION TIME WINDOW, to compute j - Deletion of the data outside this time
%interval will happen
tlim1 = 2; %in us
tlim2 = 20;%in us
% For the study case tlim1<2 us involves noise.

% Computation time window: indexes. [tlim1 tlim2] define a vector of size
% tvector:
t1= 1; % >=1, 1 corresponds to tlim1
t2= 1750; %<=length(tvector), corresponding to t2. For the ZENODO data: with tlim = [2 20]us, length(tvector) is 1800.
%This time interval will be considered for any time-series or time-integrated output.

%peak times of interest - indexes
timesofinterest= [138   334   564];     % should be consistent with your tlim1,tlim2 and the time vector used
% Other timesofinterest = [ 1 38 234 464 185 269] + 100;

%%
if computedsolution.flags.flag1==1 %j computation ON
    %----------------------------------------------------------------------------
    %----------------P R E P A R E-------T H E-------D A T A---------------------
    %----------------------------------------------------------------------------
    % 1: Load the data downloaded from ZENODO and discard from the measured data that part outside the
    % time window of interest for computation, i.e., [tlim1 tlim2].
    %----------------------------------------------------------------------------
    % Remark 1: The ZENODO data involve the time window between [0 ~30] us.
    % Remark 2: The ZENODO data are organised by distinguishing between frame positions. These frame positions are defined in agreement with the reference system presented in the
    % paper (RB is -x,-y, LB is -x,+y, LT is +x,+y and RT is +x,-y). About ZENODO data: the probe measurements within each frame position are organised indicating the relative location of each
    % probe with respect to the corresponding central one for that frame position (either H=0 or V=0),
    % i.e., +2,+4,...,+14.

    % a) Load the data
    load(motherfolder+ZENODOdata);
    % SBarquero_PPTplume_cross_section_dataset =
    %   struct with fields:
    %              t: [0 1.0000e-08 … ] %s
    %     positionRB: [1×1 struct] %A
    %     positionLB: [1×1 struct] %A
    %     positionRT: [1×1 struct] %A
    %     positionLT: [1×1 struct] %A
    database = SBarquero_PPTplume_cross_section_dataset; % more convenient
    clear SBarquero_PPTplume_cross_section_dataset % to save memory

    % b) Select the data of interest: time window between [tlim1 tlim2] %us
    tvector  =  database.t; %s
    tlims=[tlim1 tlim2]*1e-6; % tlim1,2 in us
    for t=1:length(tvector)
        if tvector(t)>tlims(1)
            tstart = t;
            break
        end
    end
    for t=tstart:length(tvector)
        if tvector(t)>=tlims(2)
            tend = t;
            break
        end
    end
    % shorten the time vector
    tvector  = tvector(tstart:tend); % s
    computedsolution.t = tvector; %this is what will be saved with the computed ion current density distribution.
    % shorten the probe data
    frame_position_list=["positionRB" "positionRT" "positionLB" "positionLT"];
    for i = 1:4
        frame_pos = frame_position_list(i);
        for j = 1:2
            if j == 1
                probe_alignment = "H";
            else
                probe_alignment = "V";
            end
            for k = 1:8
                probe_location = "plus"+string(2*(k-1));
                database.(frame_pos).(probe_alignment).means.(probe_location)=database.(frame_pos).(probe_alignment).means.(probe_location)(tstart:tend);
            end
        end
    end


    %----------------------------------------------------------------------------
    % 2: Prepare any other setup information required by the algorigthm
    %----------------------------------------------------------------------------

    % a) Setup info - PPT (1000V,6uF=3J)
    computedsolution.setupinfo.mbit = 3.9e-9; %kg //1 ug = 10^-9 kg
    % Scanner:
    requested_data.eff_coll_width   = 0.5*1e-3; %m
    Delta_x_ref                 = 2*1e-2; %m, =Delta_y_ref
    % subgrid info
    requested_data.xvector      = Delta_x_ref* ones([1 9]); %allocate

    %By default:
    computedsolution.setupinfo.gridcase = "uniform";

    if computedsolution.flags.flag1a == 1 % nonlinear system
        if computedsolution.flags.flag1a1==1 %nonuniform grid
            computedsolution.setupinfo.gridcase = "nonuniform";

            requested_data.xvector(1)   = 1.5*1e-2;
            requested_data.xvector(end) = (2 + 1.7)*1e-2;
        end
    end
    requested_data.yvector      = requested_data.xvector; % from 0 to +14
    requested_data.probelength  = sum(requested_data.xvector);

    % Plasma measurements within a subgrid. One H probe and one V probe have to be taken from the
    %corresponding neighbouring frame position.
    database.positionRB.H.means.extra_minus2 = database.positionRT.H.means.plus2 ; %minus, plus is just a reference with respect to the central probe, but has nothing to do with the reference system presented in the paper
    database.positionRB.V.means.extra_minus2 = database.positionLB.V.means.plus2 ;
    database.positionRT.H.means.extra_minus2 = database.positionRB.H.means.plus2 ;
    database.positionRT.V.means.extra_minus2 = database.positionLT.V.means.plus2 ;
    database.positionLB.H.means.extra_minus2 = database.positionLT.H.means.plus2 ;
    database.positionLB.V.means.extra_minus2 = database.positionRB.V.means.plus2 ;
    database.positionLT.H.means.extra_minus2 = database.positionLB.H.means.plus2 ;
    database.positionLT.V.means.extra_minus2 = database.positionRT.V.means.plus2 ;

    % Define the cross profiles at each t and load them into requested_data
    for i = 1:4
        frame_pos = frame_position_list(i);
        for j = 1:2
            if j == 1
                probe_alignment = "H";
            else
                probe_alignment = "V";
            end
            for k = 1:9
                if k>1
                    probe_location = "plus"+string(2*(k-2));
                else
                    probe_location = "extra_minus2";
                end
                for t= t1:t2
                    requested_data.measureddata.(frame_pos).(probe_alignment)(k,t) =  database.(frame_pos).(probe_alignment).means.(probe_location)(t); %A 
                end
            end
        end
    end
    %---------------------------------
    % For checking:
    %---------------------------------

    % figure ()
    % plot(requested_data.measureddata.positionRB.H)
    % plot(requested_data.measureddata.positionRB.V)
    % plot(requested_data.measureddata.positionLB.V)
    % plot(requested_data.measureddata.positionLB.H)
    % plot(requested_data.measureddata.positionRT.H)
    % plot(requested_data.measureddata.positionRT.V)
    % plot(requested_data.measureddata.positionLT.V)
    % plot(requested_data.measureddata.positionLT.H)
    clear database; % to save memory

    %---------------------------------

    %%
    %----------------------------------------------------------------------------
    %--C O M P U T E:---C U R R E N T---D E N S I T Y---D I S T R I B U T I O N--
    %----------------------------------------------------------------------------
    % 3: Compute for each subgrid and each t: h(x,y),v(x,y) distribution functions and J
    % a) By means of a minimizer
    % b) Analytic method - Exact solution. Only valid for uniform grid
    %----------------------------------------------------------------------------
    if computedsolution.flags.flag1a == 1

        % Algorithm information
        computedsolution.algorithmtype = "nonlinear";


        % Initial guess
        sol0(1:9)   = abs(requested_data.measureddata.positionRB.H(:,t1)*1e-1/(requested_data.eff_coll_width*1e+2 *1e+2)); %hJ, A*1e-1/cm2(perpendicular contribution assumed to be 0
        sol0(10:18) = abs(requested_data.measureddata.positionRB.V(:,t1)*1e-1/(requested_data.eff_coll_width*1e+2 *1e+2)); %vJ, A*1e-1/cm2 (perpendicular contribution assumed to be 0
        % the 1e-1 factor prepares the data for the minimizer, used later
        refJ = max (sol0); % A*1e-1/cm2
        sol0 = sol0/refJ;  %[h; v]

                %function initial guess
                sol0 = sqrt(sol0); % to prevent neg solutions
        
                A_y_overline = sum((sol0(1:9).^2).*requested_data.yvector*1e+2 );
                A_x_overline = sum((sol0(10:18).^2).*requested_data.xvector*1e+2 );
                requested_data.J = 1/(2*requested_data.eff_coll_width*1e+2)*(sum(requested_data.measureddata.positionRB.H(:,t1))*1e-1/A_y_overline  +  sum(requested_data.measureddata.positionRB.V(:,t1))*1e-1/A_x_overline ); % A*1e-1/cm2

        % Algorithm
        for i = 1:4 % for each subgrid
            frame_pos = frame_position_list(i);
            requested_data.frame_pos = frame_pos;
            for t = t1:t2
                requested_data.t = t;

                fun = @(x) nonlinearsystem(x,requested_data);
                options = optimoptions('fsolve', 'Display', 'iter', 'Diagnostics', 'on', 'MaxFunctionEvaluations', 1e5,'MaxIterations', 1e4);
                [solution, fval, exitflag, output]   = fsolve(fun, sol0, options);
                %                 fval               %diagnostics
                %                 exitflag           %diagnostics
                %                 output.iterations  %diagnostics

                h = solution(1:9).^2;
                v = solution(10:18).^2;
                computedsolution.(frame_pos).h(1:9,t)=  h; %= f(x,y,t) function in the paper
                computedsolution.(frame_pos).v(1:9,t)=  v; %= g(x,y,t) function in the paper

                % Recomputing J:
                %Expression that minimizes the error for J, as J fulfills the respective expressions for H and V (see nonlinearsystem()), but also sum(H)=sum(V)
                A_y_overline = sum(h.*requested_data.yvector*1e+2 );
                A_x_overline = sum(v.*requested_data.xvector*1e+2 );
                requested_data.J= 1/(2*requested_data.eff_coll_width*1e+2)*(sum(requested_data.measureddata.(frame_pos).H(:,t))*1e-1/A_y_overline  +  sum(requested_data.measureddata.(frame_pos).V(:,t))*1e-1/A_x_overline ); %A*1e-1/cm2
                computedsolution.(frame_pos).J(t) = requested_data.J * 1e5; %A/m2, International System, to use it outside the algorithm

                sol0  =  solution; %next initial guess

                t %time instant, for the reference of the user
            end
        end
    else
        % Algorithm information
        computedsolution.algorithmtype = "linear"; %Applies uniform cell size - The algorithm only considers xvector(2),yvector(2)
        % Algorithm
        for i = 1:4 % for each subgrid
            frame_pos = frame_position_list(i);
            requested_data.frame_pos = frame_pos;
            for t = t1:t2
                requested_data.t = t;

                [f,g,J] = linearizedsystem(requested_data); %Exact solution
                computedsolution.(frame_pos).h(1:9,t)=  f; %= f(x,y,t)  
                computedsolution.(frame_pos).v(1:9,t)=  g; %= g(x,y,t)
                computedsolution.(frame_pos).J(t) = J ; %A/m2,  
            end
        end
    end
    computedsolution.setupinfo.eff_coll_width = requested_data.eff_coll_width; %m

    %---------------------------------
    % For checking:
    %---------------------------------
    computedsolution.reconstructeddata = reconstructionfunction(computedsolution, requested_data, t1,t2,computedsolution.flags.flag1a);
    %---------------------------------

    clear requested_data; %to save memory

 


%----------------------------------------------------------------------------
% 4: Compute j(x,y,t) = h(x,y,t) v(x,y,t) J(t)
%----------------------------------------------------------------------------
% a) j within each subgrid, using the paper reference system (RB is -x,-y, LB is -x,+y, LT is +x,+y and RT is +x,-y)
% i=1,j=1 correspond to x=-14,y=-14 (RB sector of the complete grid,
% center of the outer corner cell)

%Subgrid related with the RB frame position
computedsolution.j.positionRB  = cell([1 length(tvector)]); %allocate
computedsolution.j.positionRT  = cell([1 length(tvector)]); %allocate
computedsolution.j.positionLB  = cell([1 length(tvector)]); %allocate
computedsolution.j.positionLT  = cell([1 length(tvector)]); %allocate
computedsolution.j.wholedomain = cell([1 length(tvector)]); %allocate
for t = t1:t2
    computedsolution.j.positionRB{t} = zeros([15 15]); %allocate
    computedsolution.j.positionRT{t} = zeros([15 15]); %allocate
    computedsolution.j.positionLB{t} = zeros([15 15]); %allocate
    computedsolution.j.positionLT{t} = zeros([15 15]); %allocate
    computedsolution.j.wholedomain{t}= zeros([15 15]); %allocate
    for i=1:9
        for j=1:9
            computedsolution.j.positionRB{t}(i,j) = computedsolution.positionRB.J(t)*  computedsolution.positionRB.h(10-i,t)* computedsolution.positionRB.v(10-j,t);
        end
    end
    %Subgrid related with the RT frame position
    for i=7:15
        for j=1:9
            computedsolution.j.positionRT{t}(i,j) = computedsolution.positionRT.J(t)*  computedsolution.positionRT.h(i-6,t)* computedsolution.positionRT.v(10-j,t);
        end
    end
    %Subgrid related with the LB frame position
    for i=1:9
        for j=7:15
            computedsolution.j.positionLB{t}(i,j) = computedsolution.positionLB.J(t)*  computedsolution.positionLB.h(10-i,t)* computedsolution.positionLB.v(j-6,t);
        end
    end
    %Subgrid related with the LT frame position
    for i=7:15
        for j=7:15
            computedsolution.j.positionLT{t}(i,j) = computedsolution.positionLT.J(t)*  computedsolution.positionLT.h(i-6,t)* computedsolution.positionLT.v(j-6,t);
        end
    end

    % b) j - whole grid, whole study domain
    computedsolution.j.wholedomain{t}(1:6,1:6)     = computedsolution.j.positionRB{t}(1:6,1:6);
    computedsolution.j.wholedomain{t}(10:15,1:6)   = computedsolution.j.positionRT{t}(10:15,1:6);
    computedsolution.j.wholedomain{t}(1:6,10:15)   = computedsolution.j.positionLB{t}(1:6,10:15) ;
    computedsolution.j.wholedomain{t}(10:15,10:15) = computedsolution.j.positionLT{t}(10:15,10:15);
    % central 'cross-hair':
    computedsolution.j.wholedomain{t}(7:9,1:6)     = 1/2*(computedsolution.j.positionRB{t}(7:9,1:6) + computedsolution.j.positionRT{t}(7:9,1:6) );
    computedsolution.j.wholedomain{t}(1:6,7:9)     = 1/2*(computedsolution.j.positionRB{t}(1:6,7:9) + computedsolution.j.positionLB{t}(1:6,7:9) );
    computedsolution.j.wholedomain{t}(7:9,10:15)   = 1/2*(computedsolution.j.positionLB{t}(7:9,10:15) + computedsolution.j.positionLT{t}(7:9,10:15) );
    computedsolution.j.wholedomain{t}(10:15,7:9)   = 1/2*(computedsolution.j.positionRT{t}(10:15,7:9) + computedsolution.j.positionLT{t}(10:15,7:9) );

    computedsolution.j.wholedomain{t}(7:9,7:9)     = 1/4*(computedsolution.j.positionRB{t}(7:9,7:9) + computedsolution.j.positionRT{t}(7:9,7:9) +...
    computedsolution.j.positionLB{t}(7:9,7:9) + computedsolution.j.positionLT{t}(7:9,7:9));
end

% b) Define the cell size within the complete grid for future
% calculations/postprocess
% TO PLOT:
% Reminder: i=1,j=1 correspond to x=-14,y=-14 (RB sector of the complete grid,
% center of the outer corner cell)
for i= 1:15
    computedsolution.setupinfo.coordinates.mapy_toplot(i,:) = -14:2:14;%cm
end
for j= 1:15
    computedsolution.setupinfo.coordinates.mapx_toplot(:,j) = -14:2:14;%cm
end
% TO COMPUTE. Proper cell size:
computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK      = ones([1 15]);
if computedsolution.flags.flag1a == 1 % nonlinear system
    if computedsolution.flags.flag1a1 == 1 %nonuniform grid
        computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK(1)   = 1.85; 
        computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK(end) = computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK(1);
    end
end
computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK     = computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK*Delta_x_ref; %m
computedsolution.setupinfo.coordinates.domainrow_cellsizeOK        = computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK;%m


if computedsolution.flags.flag6s==1
    save(motherfolder+computedsolution_matfile, 'computedsolution');
end

end

%%
%----------------------------------------------------------------------------
%--C O M P U T E from the C U R R E N T---D E N S I T Y---D I S T R I B U T I O N--
%----------------------------------------------------------------------------
% 6: Some computations at the 3 times of interest
%----------------------------------------------------------------------------
if computedsolution.flags.flag2==1
    % LOAD THE j DISTRIBUTION
    load(motherfolder+computedsolution_matfile)

    jmap    = computedsolution.j.wholedomain; %  A/m2
    tvector = computedsolution.t;  % s
    dt = tvector(2)-tvector(1);

    % a) Compute the cross profiles int(j(x,y,t) dx) and int(j(x,y,t) dy)
    for t=t1:t2
        for i = 1: 15
            computedsolution.j.j_projected_on_x{t}(i) = sum(jmap{t}(i,:).*computedsolution.setupinfo.coordinates.domainrow_cellsizeOK); %A/m
            j= i;
            computedsolution.j.j_projected_on_y{t}(j) = sum(jmap{t}(:,j).*computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK.');  %A/m
        end


        % b) Compute the corresponding center of mass of j
        Y = computedsolution.setupinfo.coordinates.mapx_toplot(:,1);
        X = computedsolution.setupinfo.coordinates.mapy_toplot(1,:);
        if computedsolution.flags.flag1a == 1
            if computedsolution.flags.flag1a1 == 1
                Y(1)   = -14.85; %cm
                Y(end) =  14.85;
                X(1)   = -14.85;
                X(end) =  14.85;
            end
        end
        sum_j  = sum(sum(jmap{t} ));
        center_of_mass = 0;
        for i=1:15
            for j=1:15
                center_of_mass = center_of_mass + [X(j) Y(i)].*jmap{t}(i,j);
            end
        end
        computedsolution.j.center_of_mass{t} = center_of_mass/(sum_j);
 
        %for the reference of the user
        t
        center_of_mass/(sum_j) 
    end


    %----------------------------------------------------------------------------
    % 7: Compute I_tot(t), the reconstructed time-series of the current
    % collected by the study domain and the total collected charge Qtot
    %----------------------------------------------------------------------------
    
    % a) I_tot(t)
    std_celllength = computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK(2);
    extreme_celllength = computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK(1);
    Astd = std_celllength^2;
    Aext = std_celllength*extreme_celllength;
    Acorner = extreme_celllength^2;

    for t=t1:t2

        I_tot =sum(sum(jmap{t}(2:14,2:14)))   * Astd; %A/m2 * m2

        I_tot =I_tot + sum(jmap{t}(2:14,1)) *Aext;
        I_tot =I_tot + sum(jmap{t}(2:14,15))*Aext;
        I_tot =I_tot + sum(jmap{t}(1,2:14)) *Aext;
        I_tot =I_tot + sum(jmap{t}(15,2:14))*Aext;

        I_tot =I_tot + sum(jmap{t}(1,1))  *Acorner;
        I_tot =I_tot + sum(jmap{t}(15,15))*Acorner;
        I_tot =I_tot + sum(jmap{t}(15,1)) *Acorner;
        I_tot =I_tot + sum(jmap{t}(1,15)) *Acorner;

        computedsolution.I_tot(t) = I_tot;
    end

    % b) Qtot from rho
    computedsolution.Qtot.fromI_tot = sum(computedsolution.I_tot(t1:t2)*dt); % C

    % Rough calculation to compare it with the mbit /It is assumed that the
    % grid collects all the plume current, and that the plasma is formed by
    % a single species of ions and neutrals:
    A_Carbon = 12.01; %uma
    Carbon_atom_mass = A_Carbon*1.660540199e-27 ; %kg
    Carbon_Ablatedatoms =   computedsolution.setupinfo.mbit/Carbon_atom_mass; 
    Carbon_singlechargedions_totalcharge_fullionization = Carbon_Ablatedatoms*1.6e-19; %C 
    computedsolution.m_utilization_Carbon_singlechargedions = computedsolution.Qtot.fromI_tot/Carbon_singlechargedions_totalcharge_fullionization;
    
    A_Fluorine = 17; %uma
    Fluorine_atom_mass = A_Fluorine*1.660540199e-27 ; %kg
    Fluorine_Ablatedatoms =   computedsolution.setupinfo.mbit/Fluorine_atom_mass; 
    Fluorine_singlechargedions_totalcharge_fullionization = Fluorine_Ablatedatoms*1.6e-19; %C 
    computedsolution.m_utilization_Fluorine_singlechargedions = computedsolution.Qtot.fromI_tot/Fluorine_singlechargedions_totalcharge_fullionization;

    
    %(just for the reference of the reader
    % Maximum expected avg current
    avg_dischargetime = 10e-6; %s %
    Estimated_avg_current_Carbon = Carbon_singlechargedions_totalcharge_fullionization/avg_dischargetime;                                  %5A is ok per probe
    Estimated_avg_current_Fluorine = Fluorine_singlechargedions_totalcharge_fullionization/avg_dischargetime;                                  %5A is ok per probe

 
    if computedsolution.flags.flag7s==1
        save(motherfolder+computedsolution_matfile, 'computedsolution');
    end
end

%----------------------------------------------------------------------------
% 8: Compute rho(x,y) = int (j(x,y,t) dt) within 2--20 us, its center of
% mass, cross profiles, all the charge collected by the study domain (Qtot) 
%----------------------------------------------------------------------------
if computedsolution.flags.flag3==1
    % LOAD THE j DISTRIBUTION
    load(motherfolder+computedsolution_matfile)

    jmap    = computedsolution.j.wholedomain; %  A/m2
    tvector = computedsolution.t;  % s
    dt  = tvector(2)-tvector(1);
    rho = computedsolution.j.wholedomain{1}*0; %allocate

    % a) rho(x,y)
    for i=1:15
        for j=1:15
            for t=   t1:   t2
                rho(i,j) = rho(i,j) + (jmap{t}(i,j))*dt; % jmap{t}(i,j)  in A/m2 
            end
        end
    end
    computedsolution.rho.wholedomain = rho;

    % b) Compute the corresponding center of mass of rho
    Y = computedsolution.setupinfo.coordinates.mapx_toplot(:,1);
    X = computedsolution.setupinfo.coordinates.mapy_toplot(1,:);
    if computedsolution.flags.flag1a == 1
        if computedsolution.flags.flag1a1 == 1
            Y(1)   = -14.85; %cm
            Y(end) =  14.85;
            X(1)   = -14.85;
            X(end) =  14.85;
        end
    end
    sum_rho  = sum(sum(rho ));
    center_of_mass = 0;
    for i=1:15
        for j=1:15
            center_of_mass = center_of_mass + [X(j) Y(i)].*rho(i,j);
        end
    end
    computedsolution.rho.center_of_mass{t} = center_of_mass/(sum_rho);

    %for the reference of the user
    center_of_mass/(sum_rho)

    % c) row- and column-wise profiles
    for i = 1: 15
        computedsolution.rho.rho_projected_on_x(i) = sum(rho(i,:).*computedsolution.setupinfo.coordinates.domainrow_cellsizeOK); %C/m
        j= i;
        computedsolution.rho.rho_projected_on_y(j) = sum(rho(:,j).*computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK.');  %C/m
    end

    % d) Qtot from rho
    std_celllength = computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK(2);
    extreme_celllength = computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK(1);
    Astd = std_celllength^2;
    Aext = std_celllength*extreme_celllength;
    Acorner = extreme_celllength^2;

    Qtot =sum(sum(rho(2:14,2:14))) * Astd; %C/m2 * m2

    Qtot =Qtot + sum(rho(2:14,1)) *Aext;
    Qtot =Qtot + sum(rho(2:14,15))*Aext;
    Qtot =Qtot + sum(rho(1,2:14)) *Aext;
    Qtot =Qtot + sum(rho(15,2:14))*Aext;

    Qtot =Qtot + sum(rho(15,1))   *Acorner;
    Qtot =Qtot + sum(rho(15,15))  *Acorner;
    Qtot =Qtot + sum(rho(1,1))    *Acorner;
    Qtot =Qtot + sum(rho(1,15))   *Acorner;

    computedsolution.Qtot.fromrho = Qtot; % C

    if computedsolution.flags.flag8s==1
        save(motherfolder+computedsolution_matfile, 'computedsolution');
    end
end


%----------------------------------------------------------------------------
% 9: DIVERGENCE:
% Compute instantaneous and average half-width-at-half-angle divergence angles and Qstripe/Qtot x,y profiles
% and Qstripe/Qtot along x and y
%----------------------------------------------------------------------------
% Angles are defined from the thruster discharge channel exit plane
%
if computedsolution.flags.flag4==1
    % LOAD THE j DISTRIBUTION
    load(motherfolder+computedsolution_matfile)

    jmap    = computedsolution.j.wholedomain; %  A/m2
    rho     = computedsolution.rho.wholedomain; %  C/m2
%     tvector = computedsolution.t;  % s NOT NEEDED

    % 1) Instantaneous hahw divergence angles
    half_stripe_end = [16.7  13 11 9 7 5 3 1]; %cm, with respect to (0,0), absolute values
    distance_to_exitplane = 15.3 - 4.5; % cm (distance_to_the_PTFEsurface - channellength in the PULSA APPT
    for t = t1:t2
        MAXH(t) = max(computedsolution.j.j_projected_on_x{t});% the time dependence is kept in case they are useful for internal checkings
        refH(t) = 0.5 * MAXH(t);
        for i = 1: 8 
            if computedsolution.j.j_projected_on_x{t}(i)> refH(t)
                if i>1
                    computedsolution.divergence.hahw.x_t(t)=rad2deg(atan(half_stripe_end(i-1)/distance_to_exitplane)); %deg
                else
                    computedsolution.divergence.hahw.x_t(t)=rad2deg(atan(half_stripe_end(1)/distance_to_exitplane)); %deg
                end
                break
            end
        end

        MAXV(t) =max(computedsolution.j.j_projected_on_y{t });
        refV(t) = 0.5 * MAXV(t);
        for j = 1: 8
            if computedsolution.j.j_projected_on_y{t}(j)> refV(t)
                if j>1
                    computedsolution.divergence.hahw.y_t(t)=rad2deg(atan(half_stripe_end(j-1)/distance_to_exitplane)); %deg
                else %saturated
                    computedsolution.divergence.hahw.y_t(t)=rad2deg(atan(half_stripe_end(1)/distance_to_exitplane));  %deg
                end
                break
            end
        end

    end

    % 2) Average hahw divergence angles
    MAXHrho  = max(computedsolution.rho.rho_projected_on_x);
    refHrho  = 0.5 * MAXHrho ;
    for i = 1: 8
        if computedsolution.rho.rho_projected_on_x(i)> refHrho
            if i>1
                computedsolution.divergence.hahw.x_avg =rad2deg(atan(half_stripe_end(i-1)/distance_to_exitplane)); %deg
            else
                computedsolution.divergence.hahw.x_avg =rad2deg(atan(half_stripe_end(1)/distance_to_exitplane)); %deg
            end
            break
        end
    end

    MAXVrho  = max(computedsolution.rho.rho_projected_on_y);
    refVrho  = 0.5 * MAXVrho;
    for j = 1: 8
        if computedsolution.rho.rho_projected_on_y(j)> refVrho
            if j>1
                computedsolution.divergence.hahw.y_avg =rad2deg(atan(half_stripe_end(j-1)/distance_to_exitplane));  %deg
            else %saturated
                computedsolution.divergence.hahw.y_avg =rad2deg(atan(half_stripe_end(1)/distance_to_exitplane));  %deg 
            end
            break
        end
    end


    % 3) Qstripe/Qtot along x and y
    std_celllength = computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK(2); %m
    extreme_celllength = computedsolution.setupinfo.coordinates.domaincolumn_cellsizeOK(1); %m 
    Astd = std_celllength^2; %m2
    Aext = std_celllength*extreme_celllength;%m2
    Acorner = extreme_celllength^2;%m2

    computedsolution.setupinfo.coordinates.domaincolumn_cellcentercoordOK       = [0 2 4 6 8 10 12 15.85] * 0.01; %m
    computedsolution.setupinfo.coordinates.domainrow_cellcentercoordOK          = computedsolution.setupinfo.coordinates.domaincolumn_cellcentercoordOK ; %m
    
    computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x            = computedsolution.setupinfo.coordinates.domaincolumn_cellcentercoordOK;
    computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x(2:end-1)   = computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x(2:end-1) -0.01;%m
    computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x(end)       = 13*1e-2;%m
    computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x(length(computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x) + 1) = 16.7*1e-2;%m
    computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x = computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x; %m
    computedsolution.setupinfo.coordinates.stairs_spacecoordinates_y = computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x; %m

    Qstripe_x = zeros([1 8]);
    Qstripe_y = zeros([1 8]);
    for i = 8: 14
        j=i;
        %x
        Qstripe = sum(sum(rho(8-(i-8):8+(i-8),2:14)))*Astd ;
        Qstripe = Qstripe + sum(sum(rho(8-(i-8):8+(i-8),1)))*Aext  ;
        Qstripe_x(i-7) = Qstripe + sum(sum(rho(8-(i-8):8+(i-8),15)))*Aext  ;
 
        %y
        Qstripe = sum(sum(rho(2:14, 8-(j-8):8+(j-8))))*Astd ;
        Qstripe = Qstripe + sum(sum(rho(1,8-(i-8):8+(i-8))))*Aext  ;
        Qstripe_y(j-7) = Qstripe + sum(sum(rho(15,8-(i-8):8+(i-8))))*Aext  ;
 
    end
    Qstripe_x(8) = computedsolution.Qtot.fromrho; %C
    Qstripe_y(8) =  Qstripe_x(8);

    computedsolution.divergence.Qratio_x = (Qstripe_x/computedsolution.Qtot.fromrho)*100; % \%
    computedsolution.divergence.Qratio_y = (Qstripe_y/computedsolution.Qtot.fromrho)*100; % \%

    if computedsolution.flags.flag9s==1
        save(motherfolder+computedsolution_matfile, 'computedsolution');
    end
end



%%
%----------------------------------------------------------------------------
%---------------P O S T P R O C E S S----------------------------------------
%----------------------------------------------------------------------------


if computedsolution.flags.flag10 ==1

    % LOAD THE j DISTRIBUTION
    load(motherfolder+ZENODOdata)
    database = SBarquero_PPTplume_cross_section_dataset; % to be more convenient
    clear SBarquero_PPTplume_cross_section_dataset % to save memory
    
    load(motherfolder+computedsolution_matfile)

    %Distribution maps
    jmap    = computedsolution.j.wholedomain; %  A/m2
    rho     = computedsolution.rho.wholedomain; %  C/m2
    tvector = computedsolution.t;  % s NOT NEEDED

    %Grid coordinates
    mapy_toplot = computedsolution.setupinfo.coordinates.mapy_toplot;
    mapx_toplot = computedsolution.setupinfo.coordinates.mapx_toplot;

    %----------------------------------------------------------------------------
    % 1. Plot j(x,y) at the 3 time of interest
    %----------------------------------------------------------------------------
    for t= timesofinterest
        max_j        = max(max(computedsolution.j.wholedomain{t}));
        MAX(t)       = max (max_j) ;
    end
    MAX = max(MAX);

    for t= timesofinterest

        figure()
        map=bar3(rot90((jmap{t}.') ), 1)

        %--------------------------------------------------------------------------
        %Format  
        %--------------------------------------------------------------------------
        title('t='+string( tvector(t)*1e6)+' us, index='+string(t))

        %     map.EdgeColor = 'none';

        colormap(flipud(turbo))
        A         = colorbar();

        A.Limits  = [0 MAX];

        A.Label.String      = '$j$,  Am$^{-2}$';
        A.Label.FontSize    = 14;
        A.Label.Interpreter = 'latex';
        %     A.Label.Rotation = 270;

        xlabel('$y$, cm', 'Interpreter','latex', 'fontsize', 14)
        ylabel('$x$, cm', 'Interpreter','latex', 'fontsize', 14)
        zlabel ('$j$,  Am$^{-2}$', 'Interpreter','latex', 'fontsize', 14,'Rotation',270)
        set(gca, 'XMinorTick','on','YMinorTick','on','fontsize', 14)

        view(2)

        caxis([0 MAX])
        zlim([0 MAX])
        ylim([0.5  15.5])
        xlim([0.5  15.5])

        yticks([1:1:+15])
        xticklabels([-14:+2:+14])
        yticklabels([14:-2:-14])
        xticks([1:1:+15])

        axis square

        %Interpolation
        for map_k  = 1:length(map)
            zdata  = map(map_k).ZData;
            map(map_k).CData     = zdata;
            map(map_k).FaceColor = 'interp';
        end
        %--------------------------------------------------------------------------
        if computedsolution.flags.flag11s == 1
            savefig(gcf, motherfolder+'figure1_jmap_tindex'+string(t)+'.fig' )
        end

        %----------------------------------------------------------------------------
        % 2. Plot jprofiles at the 3 time of interest
        %----------------------------------------------------------------------------
        figure()

        subplot(1,2,1)
        stairs(-14:2:+16,[computedsolution.j.j_projected_on_y{t} 0], 'r', 'LineWidth', 2)
        hold on
        ylabel('j line integral along x, A/m')
        xlabel('y, cm')

        subplot(1,2,2)
        stairs([0 computedsolution.j.j_projected_on_x{t}],-14:2:+16, 'r', 'LineWidth', 2)
        xlabel('j line integral along y, A/m')
        xlabel('x, cm')

        sgtitle('t='+string( tvector(t)*1e6)+' us, index='+string(t))
       
        if computedsolution.flags.flag11s == 1
            savefig(gcf, motherfolder+'figure2_jprofiles_tindex'+string(t)+'.fig' )
        end
    end
    %----------------------------------------------------------------------------
    % 3. Plot rho(x,y)
    %----------------------------------------------------------------------------

    figure()
    map=bar3(rot90((rho.') ), 1)

    %Format
    %--------------------------------------------------------------------------

    %     map.EdgeColor = 'none';

    colormap(flipud(gray))
    A         = colorbar();

    max_rho     = max(max(computedsolution.rho.wholedomain));
    MAX       = max (max_rho) ;
    A.Limits  = [0 MAX];

    A.Label.String      = '$rho$,  Cm$^{-2}$';
    A.Label.FontSize    = 14;
    A.Label.Interpreter = 'latex';
    %     A.Label.Rotation = 270;

    xlabel('$y$, cm', 'Interpreter','latex', 'fontsize', 14)
    ylabel('$x$, cm', 'Interpreter','latex', 'fontsize', 14)
    zlabel ('$rho$,  Cm$^{-2}$', 'Interpreter','latex', 'fontsize', 14,'Rotation',270)
    set(gca, 'XMinorTick','on','YMinorTick','on','fontsize', 14)

    view(2)

    caxis([0 MAX])
    zlim([0 MAX])
    ylim([0.5  15.5])
    xlim([0.5  15.5])

    yticks([1:1:+15])
    xticklabels([-14:+2:+14])
    yticklabels([14:-2:-14])
    xticks([1:1:+15])

    axis square

    %Interpolation
    for map_k  = 1:length(map)
        zdata  = map(map_k).ZData;
        map(map_k).CData     = zdata;
        map(map_k).FaceColor = 'interp';
    end

    %--------------------------------------------------------------------------
    if computedsolution.flags.flag11s == 1
        savefig(gcf, motherfolder+'figure3_rhomap.fig' )
    end

    %----------------------------------------------------------------------------
    % 4. Plot rhoprofiles
    %----------------------------------------------------------------------------
    figure()

    subplot(1,2,1)
    stairs(-14:2:+16,[computedsolution.rho.rho_projected_on_y 0], 'r', 'LineWidth', 2)
    hold on
    ylabel('rho line integral along x, C/m')
    xlabel('y, cm')

    subplot(1,2,2)
    stairs([0 computedsolution.rho.rho_projected_on_x],-14:2:+16, 'r', 'LineWidth', 2)
    xlabel('rho line integral along y, C/m')
    xlabel('x, cm')

    if computedsolution.flags.flag11s == 1
        savefig(gcf, motherfolder+'figure4_rhoprofiles.fig' )
    end

    %----------------------------------------------------------------------------
    % 5. Plot I_tot(t)
    %----------------------------------------------------------------------------
    figure()

    plot(tvector(t1:t2)*1e6, computedsolution.I_tot(t1:t2),'Color',  [0 0 0], 'LineWidth',2)


    xlabel('$t$, $\mu$s', 'interpreter', 'latex', 'FontSize',14)
    ylabel('$I_{tot}$, A', 'interpreter', 'latex', 'FontSize',14)
    set(gca, 'XMinorTick','on','YMinorTick','on','fontsize', 14)
    xlim([2 20])

    grid on
    axis square


    if computedsolution.flags.flag11s == 1
        savefig(gcf, motherfolder+'figure5_I_tot.fig' )
    end

    %----------------------------------------------------------------------------
    % 6. Qratio
    %----------------------------------------------------------------------------

    stairs_spacecoordinates_x=computedsolution.setupinfo.coordinates.stairs_spacecoordinates_x*1e+2; %cm
    stairs_spacecoordinates_y=stairs_spacecoordinates_x; %cm

    figure()

    stairs (stairs_spacecoordinates_x,[computedsolution.divergence.Qratio_x computedsolution.divergence.Qratio_x(end)], 'k', 'linewidth', 2)
    hold on
    stairs (stairs_spacecoordinates_y,[computedsolution.divergence.Qratio_y computedsolution.divergence.Qratio_y(end)], 'r', 'linewidth', 2)
    ylim([0 100])
    xlim([0 16.7])
    axis square
    ylabel('$Q_{stripe}/Q_{tot}$, %','interpreter', 'latex','FontSize',14)
    xlabel ('$x$ (black) or $y$ (red), cm','interpreter', 'latex','FontSize',14 )


    if computedsolution.flags.flag11s == 1
        savefig(gcf, motherfolder+'figure6_divergence_Qratio.fig' )
    end

    %----------------------------------------------------------------------------
    % 7. hahw divergence angles: instantaneous and average
    %----------------------------------------------------------------------------
    figure()

    plot(tvector(t1:t2)*1e6,computedsolution.divergence.hahw.x_t, 'k', 'LineWidth',2)
    hold on
    plot(tvector(t1:t2)*1e6,computedsolution.divergence.hahw.y_t, 'r', 'LineWidth',2)
    hold on
    xlim([2 20])

    plot(tvector(t1:t2)*1e6,computedsolution.divergence.hahw.x_avg * ones ([ t2 1]), 'k:', 'LineWidth',2)
    hold on
    plot(tvector(t1:t2)*1e6,computedsolution.divergence.hahw.y_avg * ones ([ t2 1]), 'r:', 'LineWidth',2)
    hold on

    ylabel('$\theta$, $^\circ$',  'interpreter', 'latex', 'FontSize',14)
    xlabel('$t$, $\mu$s',  'interpreter', 'latex', 'FontSize', 14)

    legend({'theta_x(t)' 'theta_y(t)' 'theta_{x,avg}' 'theta_{y,avg}'})

    xlim([2 15])

    axis square

    if computedsolution.flags.flag11s == 1
        savefig(gcf, motherfolder+'figure7_divergence_hahw.fig' )
    end


    if computedsolution.flags.flag00 == 1
        %----------------------------------------------------------------------------
        % Additional postprocess for debugging -
        %----------------------------------------------------------------------------
        % 1)Reconstruction of measured data
        %----------------------------------------------------------------------------
        frame_position_list=["positionRB" "positionRT" "positionLB" "positionLT"];
        alignment_list = ["H" "V"];
        for k = 1:4
            frame_pos = frame_position_list(k);

            figure ()

            for a=1:2
                subplot(1,2,a)

                alignment = alignment_list (a);
                for l = [0 2 4 6 8 10 12 14]
                    probelocation = "plus"+string(l);

                    plot(tvector(t1:t2)*1e6,computedsolution.reconstructeddata.(alignment).(frame_pos).(probelocation))
                    hold on

                    plot(database.t*1e6,database.(frame_pos).(alignment).means.(probelocation), 'r:')
                    hold on
                end
                xlim([2 15])
                ylim([0 9])

                xlabel('t, us')
                ylabel(alignment+', A')
            end
            sgtitle(frame_pos)


            if computedsolution.flags.flag11s == 1
                if computedsolution.flags.flag1a == 1
                    % b) Nonlinear system
                    meshtype = computedsolution.setupinfo.gridcase;
                    savefig(gcf, motherfolder+'figure8_'+frame_pos+'_checks_reconstruction_nonlinearized_'+meshtype+'.fig' )
                else
                    % a) Linearized system
                    savefig(gcf, motherfolder+'figure8_'+frame_pos+'_checks_reconstruction_linearized.fig' )
                end
            end
        end

        %----------------------------------------------------------------------------
        % 1) J(t) for each frame position
        %----------------------------------------------------------------------------
        figure()
        for k = 1:4
            frame_pos = frame_position_list(k);
            plot(tvector(t1:t2)*1e6,computedsolution.(frame_pos).J)
            hold on
        end
        xlim([2 15])
        xlabel('t, us')
        ylabel('J, A/m2')
        legend(frame_position_list)

        if computedsolution.flags.flag11s == 1
            if computedsolution.flags.flag1a == 1
                % b) Nonlinear system
                meshtype = computedsolution.setupinfo.gridcase;
                savefig(gcf, motherfolder+'figure9_checks_J_nonlinearized_'+meshtype+'.fig' )
            else
                % a) Linearized system
                savefig(gcf, motherfolder+'figure9_checks_J_linearized.fig' )
            end
        end
    end
end



%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
 


%% Algorithm functions

function fun0 = nonlinearsystem(x,requested_data)

frame_pos = requested_data.frame_pos;
t = requested_data.t;

% Initial guess
h = x(1:9).^2;
v = x(10:18).^2;
J = requested_data.J; %A/cm^2 *1e-1

data.H = requested_data.measureddata.(frame_pos).H * 1e-1; %A* 0.1
data.V = requested_data.measureddata.(frame_pos).V * 1e-1; %A* 0.1

fun0 = x*0; %allocate

% Appropiate units for the minimizer
eff_coll_width = requested_data.eff_coll_width *1e+2; %cm
yvector    = requested_data.yvector *1e+2; %cm
xvector    = requested_data.xvector *1e+2; %cm


% Equations 5 of section III of the paper: CROSS-SECTION RECONSTRUCTION
% ALGORITHM
% Horizontal:
for probe = 1: 9
    perpendicular_contribution = 0; %initialize
    for perpendicular_probe = 1: 9
        perpendicular_contribution = perpendicular_contribution + v(perpendicular_probe)* yvector(perpendicular_probe);
    end
    fun0(probe) = J *  eff_coll_width * h(probe)*perpendicular_contribution - data.H(probe,t);
end
% Vertical:
for probe = 1: 9
    perpendicular_contribution = 0; %initialize
    for perpendicular_probe = 1: 9
        perpendicular_contribution = perpendicular_contribution + h(perpendicular_probe)* xvector(perpendicular_probe);
    end
    fun0(probe+9) = J *  eff_coll_width * v(probe)*perpendicular_contribution - data.V(probe,t);
end

% fun0(19)= 1/(2*eff_coll_width)*(sum(data.H(:,t))/A_y_overline  +  sum(data.V(:,t))/A_x_overline ) - J;
fun0(19)= sum(h) - 1;
fun0(20)= sum(v) - 1;

% fun0(19)
% sum(fun0)

end

function [f,g,J] = linearizedsystem(requested_data)


% --------------------------------------------------------------------------

% Normalization applied to h,v:
% sum(h.*xvector)=Lx;
% sum(v.*xvector)=Ly;
% where Lx, Ly are the respective probe lengths

%Grid: 15x15cm2, and it is assumed that no current is collected by the
%probe extreme extending beyond the 15x15cm2 study domain.

% --------------------------------------------------------------------------


%Info

frame_pos = requested_data.frame_pos;
t = requested_data.t;

Deltax   = requested_data.xvector(2) ; % m
Deltay   = requested_data.yvector(2) ; % m

Imeas_H   = requested_data.measureddata.(frame_pos).H(:,t) ; %A
Imeas_V   = requested_data.measureddata.(frame_pos).V(:,t) ; %A

w       = requested_data.eff_coll_width ; %  m

Lx =  Deltax*9; % m
Ly =  Deltay*9; % m
    
% Core - exact solution

J = 1/(2*Lx*Ly*w)*((sum(Imeas_H)*Deltax)  + (sum(Imeas_V)*Deltay) ); %Expression that minimizes the error for J, as J fulfills the respective expressions for H and V (below), but also sum(H)=sum(V)
 
f = zeros([1 9]); %allocate
g = f; % allocate
for i=1:9
    if Imeas_H(i) <0
        Imeas_H(i) = 0;
    end
    f(i) = Imeas_H(i)/ (J*w*Ly); %=h
end
for j=1:9
    if Imeas_V(j) <0
        Imeas_V(j) = 0;
    end
    g(j) = Imeas_V(j)/ (J*w*Lx); %=v
end

    
end

function reconstructeddata = reconstructionfunction(computedsolution, requested_data, t1,t2, flag1a)
frame_position_list=["positionRB" "positionRT" "positionLB" "positionLT"];

for k = 1:4
    frame_pos = frame_position_list(k);
    for t=t1:t2

        cross_contribution_v = 0;
        cross_contribution_h = 0;

        if  flag1a ==1

            for i=1:9
                cross_contribution_v = cross_contribution_v + computedsolution.(frame_pos).v(i,t)*requested_data.yvector(i) ;
                cross_contribution_h = cross_contribution_h + computedsolution.(frame_pos).h(i,t)*requested_data.xvector(i) ;
            end

        else
            Lx =  requested_data.xvector(2)*9; % m
            Ly =  requested_data.yvector(2)*9; % m
            cross_contribution_v = Ly;
            cross_contribution_h = Lx;
        end

        count=2;
        for l = [0 2 4 6 8 10 12 14]
            probelocation = "plus"+string(l);
            reconstructeddata.H.(frame_pos).(probelocation)(t) = computedsolution.(frame_pos).J(t)*computedsolution.(frame_pos).h(count,t)*requested_data.eff_coll_width*(cross_contribution_v);
            reconstructeddata.V.(frame_pos).(probelocation)(t) = computedsolution.(frame_pos).J(t)*computedsolution.(frame_pos).v(count,t)*requested_data.eff_coll_width*(cross_contribution_h);
            count = count+1;
        end

    end
end
end



 
