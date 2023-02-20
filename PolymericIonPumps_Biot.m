% Created 02/15/2022 - Code for JMS Polymeric Ion Pumps Manuscript -
% version: Matlab 2022b

% Creator - Jonathan Aubuchon Ouimet
% PI - William A. Phillip
% Contributing Author - Alexander Dowling


clear
clc
clf

% We are interested in solving for the pseude steady state unloading and
% loading concentration profiles. In order to do so, we use a numerical 
% solving technique that given initial conditions (i.e., fully loaded or empty)
% oscillates between the loading and unloading stages in order
% to converge to PSS. 

% Material Properties are defined...
Bi = 0;
Cx = 100;
% System Operating Conditions are defined...
tul = 0.5;
tload_cycle = 0.4448;
% The number is discrete elements that the membrane is divided into is
% defined...
Discritization = 10001;
    
% We define the eigen values within vector 'sig'. They are dependent only on the biot number. 
sig = ones([20 1]); % A vector of ones is established to store 20 eigen values.
for i = 0:19  % The 'for' loop establishes 20 eigen values. 
            syms s % The variable s is used to write the eigen value equation. 
            sig(i+1,1) = vpasolve(-s*cot(s)-Bi==0,s,[-0.0001+(i)*-pi (i+1)*-pi+0.0001]); % 'vpasolve' is used to find the eigen values. The eigen values must be negative and the search range is changed to capture the subsequent eigen value. 
end
     

% Initialize the system and calculate integration constants.      
xi = linspace(0,1,Discritization)'; % Discretize the membrane.
for Cycle = 1:1000
    if Cycle == 1 % For the first cycle, we assume fully loaded (i.e., Cm = Cx) or unloaded profile (Cm = 0).
        Cs = 0;
    elseif Cycle > 1 % For every subsequent cycle, we use the previous loading profile. 
        Cs = LoadingProfile(:,2);
    end 
     
    % The integration constants for this cycle are calaculated.  
    IntConstant = ones([19 1]); % The integration constants are initalized within a vector of ones.
    for Int = 1:20
        % Compute values for the Eq. #S20
        % Note that ctransient = c-css. 

        Eq_top = (Cs-(Bi/(Bi+1).*(1-xi))).*(cos(sig(Int).*xi)+Bi/sig(Int).*sin(sig(Int).*xi));
        Int_top = cumtrapz(xi, Eq_top); % Integrate equation with "cumtrapz"
        Eq_bot = (cos(sig(Int)*xi)+Bi/sig(Int)*sin(sig(Int)*xi)).^2;
        Int_bot = cumtrapz(xi, Eq_bot); % Integrate equation with "cumtrapz"
        IntConstant(Int,1) = Int_top(Discritization)/Int_bot(Discritization);
    end

    % We apply the unloading stage equation to each finite element. 
    UnloadingProfile = ones([Discritization 2]); % The 1001 segments are initalized within a vector of 'ones'
    for Loop_Position = 1:Discritization
        xi_UnloadingLoop = xi(Loop_Position);
        ex = exp(-(sig.^2).*tul);
        TransConc = IntConstant.*ex.*(cos(sig*xi_UnloadingLoop)+Bi./sig.*sin(sig*xi_UnloadingLoop));
        SSConc = Bi/(Bi+1)*(1-xi_UnloadingLoop); 
        Conc = SSConc + sum(TransConc);

        UnloadingProfile(Loop_Position,1) = xi_UnloadingLoop;
        UnloadingProfile(Loop_Position,2) = Conc;
    end

    % Several key metrics to determine whether or not the system is at PSS
    % are determined. 
    ExitFlux = (Bi./(Bi+1)).*tul+sum((1-ex).*IntConstant./(sig.^2).*(sig.*sin(sig)-Bi.*cos(sig))); % The flux into the recieving solution.
    EntranceFlux = (Bi./(Bi+1)).*tul - sum(Bi*IntConstant./sig.^2.*(1-ex)); % The flux into the feed.
    TotalFluxOut = ExitFlux + abs(EntranceFlux); % The total flux out of the sorbent matrix. 
    Performance = ExitFlux/(tul+tload_cycle); % The performance for this cycle.
    

    % Now, the penetration front is solved for. 
    TloadInt = (Cx-UnloadingProfile(:,2)).*UnloadingProfile(:,1); % Define the equation that describes the penetraiton front.
    io = interp1(cumtrapz(xi,TloadInt),UnloadingProfile(:,1),tload_cycle); % Integrate previous equation and determine what the penetration front is. 

    % We now need to define the loading profile. This will be a
    % step wise function.
    LoadingProfile = ones([Discritization 2]); % The 1001 segments are initalized within a vector of 'ones'
    for LoadingLoop_Position = 1:Discritization
        xi_LoadingLoop = xi(LoadingLoop_Position);
        
        LoadingProfile(LoadingLoop_Position,1) = xi_LoadingLoop;
        if xi_LoadingLoop <= io
            LoadingProfile(LoadingLoop_Position,2) = Cx;
        elseif xi_LoadingLoop > io
            LoadingProfile(LoadingLoop_Position,2) = UnloadingProfile(LoadingLoop_Position,2);
        end
    end
    % The total flux into the membrane is a key metric to determine whether of not the system is at PSS.  
    TotalFluxIn = interp1(UnloadingProfile(:,1),cumtrapz(xi,Cx-UnloadingProfile(:,2)),io);

    % The cycle is finished and we move on to the next cycle. 
end

% The final unloading stage concentration profile is graphed. 
plot(UnloadingProfile(:,1),UnloadingProfile(:,2))