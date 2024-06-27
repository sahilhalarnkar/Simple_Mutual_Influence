[sel, enh,C,theta,theta_noc,kf,kr] = pareto_enhsel(0.01, 1, 0.01);
figure
plot(sel,enh)

% Ensure that all points have converged completely.
% Filter out points that haven't converged provided that there are enough neighbouring points to ensure an uninterrupted continuous boundary.

%%

function [enh, C, theta,theta_noc, kf, kr] = opt_enhance(eps)
    
    % Creating optimization problem seeking to maximize
    rxn_selvsenh = optimproblem('ObjectiveSense', 'maximize');
    
    % Create optimization options
    options = optimoptions(@fmincon,'TolF',1e-8);

    % Initial guesses for rate constants 
    kf_guess = [1, 2, 1, 1.2];
    kr_guess = [0.1, 0.1];
    % Define variables and auxiliary coverages
    % Define bounds to span a range of acceptable rate parameters
   
    kf = optimvar('kf', 4,1, 'LowerBound', 0.01, 'UpperBound', 10);
    kr = optimvar('kr', 2,1, 'LowerBound', 0.01, 'UpperBound', 10);
    C = optimvar('C', 3,1, 'LowerBound', 0.01, 'UpperBound', 1);
    theta = optimvar('theta', 2,1, 'LowerBound', 0.01, 'UpperBound', 1);
    theta_noc = optimvar('theta_noc', 2,1, 'LowerBound', 0.01, 'UpperBound', 1);

    % Create optimization problem
    % The objective function is the enhancement factor EF
    rxn_selvsenh.Objective=(theta(1)*theta(2))/(theta_noc(1)*theta_noc(2));

    % Constraints 
    % Summation of fractional coverages of R1, R2 and H should be less than or equal to one
    % The problem is optimized at a constant selectivity using an epsilon
    % constraint
    coverage = sum(theta) <= 1;
    coverage_noc = sum(theta_noc) <= 1;
    select_eps = kf(3)*theta(2) - eps*(kf(3)*theta(2)+kf(4)*C(3))==0;
    % Nonlinear d/dt equations
    dBdt = kf(1)*C(2)*(1-sum(theta))-kr(1)*theta(2)-kf(3)*theta(2)*theta(1) == 0;
    dAdt = kf(2)*C(1)*(1-sum(theta))-kr(2)*theta(1)-kf(3)*theta(2)*theta(1)-kf(4)*C(3)*theta(1) == 0;
    dBdt_noc = kf(1)*C(2)*(1-sum(theta_noc))-kr(1)*theta_noc(2)-kf(3)*theta_noc(2)*theta_noc(1) == 0;
    dAdt_noc = kf(2)*C(1)*(1-sum(theta_noc))-kr(2)*theta_noc(1)-kf(3)*theta_noc(2)*theta_noc(1) == 0;


    % Add constraints to the problem
    rxn_selvsenh.Constraints.coverage = coverage;
    rxn_selvsenh.Constraints.coverage_noc = coverage_noc;
    rxn_selvsenh.Constraints.dBdt = dBdt;
 
    rxn_selvsenh.Constraints.dAdt = dAdt;
    rxn_selvsenh.Constraints.dBdt_noc = dBdt_noc;
    rxn_selvsenh.Constraints.dAdt_noc = dAdt_noc;
    rxn_selvsenh.Constraints.select_eps = select_eps;

    % Solve the problem using the following initial values
    x0.kf=[10, 10, 10, 10];
    x0.kr=[10, 10];
    x0.C=[1,1,1];
    x0.theta=[1,1];
    x0.theta_noc=[1,1];

    [sol, fval, exitflag] = solve(rxn_selvsenh,x0,'Options',options);

    % Retrieve solution
    enh = fval;
    C = sol.C;
    theta = sol.theta;
    kf = sol.kf;
    kr = sol.kr;
    theta_noc=sol.theta_noc;
end

function [sel, enh,C,theta,theta_noc,kf,kr] = pareto_enhsel(selmin, selmax, delta)
    % Function pareto_enhsel
        %
        % Parameters: 
        % 
        % selmin - Lower bound of selectivity
        % selmax - Upper bound of selectivity
        % delta  - Step Size
        %
        % Returns the following variables for the optimized values of EF:
        %
        % sel: Selectivity
        % enh: Enhancement Factor
        % C : Concentrations of R1, R2 and H 
        % theta: Coverages for the R1, R2 and H system
        % theta_noc: Coverages for the R1 and H system
        % kf: Forward reaction rate constants 
        % kr: Reverse reaction rate constants 
    sel = selmin:delta:selmax;
    enh = zeros(size(sel));
    C = zeros(numel(sel),3);
    theta = zeros(numel(sel),2);
    theta_noc = zeros(numel(sel),2);
    kf = zeros(numel(sel),4);
    kr = zeros(numel(sel),2);

    % Loop through different selectivities to generate Pareto frontier
    for i = 1:numel(sel)
        [enh(i),C(i,:),theta(i,:),theta_noc(i,:),kf(i,:),kr(i,:)] = opt_enhance(sel(i));
    end
   
end



