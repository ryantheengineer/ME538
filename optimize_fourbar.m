function [xopt, fopt, exitflag, output] = optimize_fourbar()


    % ------------Starting point and bounds------------
    %var= r1  r2  r3  r4  r5  dy  a3  b3  K2 K3 theta20 theta30  %design variables
    x0 = [10, 3]; %starting point
    ub = [50, 5]; %upper bound
    lb = [5,  1]; %lower bound

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        
        %design variables
        r1 = x(1);          % Link 1 length
        r2 = x(2);          % Link 2 length
        r3 = x(3);          % Link 3 length
        r4 = x(4);          % Link 4 length
        r5 = x(5);          % Link 5 length
        dy = x(6);          % Offset between link 1 and slider of link 5
        a3 = x(7);          % Distance to actuation point along link 3
        b3 = x(8);          % Distance to actuation point perp. to link 3
        K2 = x(9);          % Link 2 spring constant
        K3 = x(10);         % Link 3 spring constant
        theta20 = x(11);    % Resting angle of link 2
        theta30 = x(12);    % Resting angle of link 3
        
        
        %other analysis variables
        theta2 = linspace(0,pi,200);
        closeangle = 0;
        openangle = 80; % Degrees (for simplicity)
        openangle = openangle*pi/180;
        
        % Binary indicators of which angle data is within the desired
        % actuation range. Use as theta2(goodrange) or V(goodrange) later
        % on. This will allow force searching within only the expected
        % angles if desired.
        goodrange = (theta2<openangle);
        
        
        
        %% Analysis Functions
        % Potential energy curve
        V = get_potential_energy(r1,r2,r3,r4,theta2,theta20,theta30,K2,K3);
        TF = islocalmin(real(V));
        nmins = nnz(TF);
        if nmins ~= 2
            I = zeros(2,1);
            [min1,I(1)] = min(real(V));
            I(2) = find(TF);
            min2 = real(V(I(2)));
        else
            I = find(TF);
            min1 = real(V(I(1)));
            min2 = real(V(I(2))); 
        end
        
        % Calculate the input force at the 
        
%         
%         weight = rho*2*pi*d*t*sqrt((B/2)^2+H^2); %lbs
%         stress = (P*sqrt((B/2)^2+H^2))/(2*t*pi*d*H); %ksi
%         bstress = (pi^2*E*(d^2+t^2))/(8*((B/2)^2+H^2)); %ksi
%         deflection = P*((B/2)^2+H^2)^(3/2)/(2*t*pi*d*H^2*E); %in

        %objective function
        f = weight; %minimize weight
        
        %inequality constraints (c<=0)
        c = zeros(3,1);         % create column vector
        c(1) = stress-100;      %stress <= 100ksi
        c(2) = stress-bstress;  %stress-bstress <= 0
        c(3) = deflection-.25;  %deflection <= 0.25
        
        %equality constraints (ceq=0)
        ceq = [];

        
    end


    % ------------Call fmincon------------
    options = optimoptions(@fmincon,'display','iter-detailed','Diagnostics','on');
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
    xopt %design variables at the minimum
    fopt %objective function value at the minumum  fopt = f(xopt)
    

    % ------------Separate obj/con (do not change)------------
    function [f] = obj(x) 
        [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x) 
        [~, c, ceq] = objcon(x);
    end
end