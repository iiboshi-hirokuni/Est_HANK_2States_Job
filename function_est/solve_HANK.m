function [G1, impact,inv_state_red,from_spline, eu, rrho] = solve_HANK(I,J,n_v,n_g,n_p,n_shocks,para)

warning('off','AutoDiff:maxmin')

rrho=0;

%%
% Set options for this example run
ReduceDistribution = 1;  % 1 for state space reduction 0 for not
reduceV = 1;             % 1 for value function reduction 0 for not
ReduceDist_hor = 20;     % Dimensionality of the Krylov subspace
DisplayLev = 1;          % Determines verbosity of steady state calculation
check_consistency = 0;   % Runs Step 6: Internal consistency check


%% Step 0: Set parameters
% The script sets up parameters relevant for the model
%    Economic parameters, approximation parameters, and grid parameters are
%    defined in the script.
set_parameters;

parameters_update;

% n_v = I*J + 1;    % number of jump variables (value function + inflation)
% n_g = I*J + 2 ;   % number of endogeneous state variables (distribution + monetary + Fiscal policy)
% n_p = 6;          % number of static relations: bond-market clearing, labor market clearing, consumption, output, total assets
% n_shocks = 3;     % only monetary policy shock is considered
nEErrors = n_v;
nVars = n_v + n_g + n_p;

%% Step 1: Solve for the Steady State
% Any methods can be used to solved for the steady state. In particular, example
%    codes can be found at \<http://www.princeton.edu/~moll/HACTproject.html\>.
fprintf('Computing steady state...\n');
% t0 = tic;
compute_steady_state;
% fprintf('Time to compute steady state: %2.4f seconds\n\n\n',toc(t0));


%% Step 2: Linearize Model Equations
% For computing derivatives, the codes written for solving for the
%    steady-state can be used almost verbatim using automatic
%    differentiation toolbox as long as only the functions supported by
%    automatic differentation are used. For list of supported functions and
%    documentation of relevant syntax check
%    <https://github.com/sehyoun/MATLABAutoDiff>. Example usage/syntax of
%    automatic differentiation can be found at
%    <https://sehyoun.com/EXAMPLE_AutoDiff_Syntax.html>
% fprintf('Taking derivatives of equilibrium conditions...\n');
%  t0 = tic;

% Prepare automatic differentiation
vars = zeros(nVars + nVars + nEErrors + n_shocks,1);
vars = myAD(vars);

% Evaluate derivatives
equilibrium_conditions;

% Extract derivative values
derivs = getderivs(v_residual);


%% Step 3: Solve out Static Constraints or Reduce the Model
% Extract derivatives
g1 = -derivs(:,1:nVars);
g0 = derivs(:,nVars+1:2*nVars);
pi = -derivs(:,2*nVars+1:2*nVars+nEErrors);
psi = -derivs(:,2*nVars+nEErrors+1:2*nVars+nEErrors+n_shocks);
constant = sparse(nVars,1);

% State Variables
if ReduceDistribution	== 1
    % Reduce model
% 	fprintf('Reducing distribution ...\n');
    [state_red,inv_state_red,n_g_red] = krylov_reduction(g0,g1,n_v,n_g,ReduceDist_hor);
    [g1,psi,pi,constant,g0] = change_basis(state_red,inv_state_red,g1,psi,pi,constant,g0);
else
    % Solve out static constraints
    fprintf('Solving Out Static Constraints ...\n');
    [state_red,inv_state_red,g0,g1,constant,pi,psi] = clean_G0_sparse(g0,g1,constant,pi,psi);
     n_g_red = n_g;
end

% Jump Variables
if reduceV == 1
    % Reduce dimensionality of value function using splines
    n_knots = 12;
    c_power = 1;
    x = a';
    n_post = size(z,2);
    n_prior = 1;

    % Create knot points for spline (the knot points are not uniformly spaced)
    knots = linspace(amin,amax,n_knots-1)';
    knots = (amax-amin)/(2^c_power-1)*((knots-amin)/(amax-amin)+1).^c_power+ amin-(amax-amin)/(2^c_power-1);

    % Function calls to create basis reduction
    [from_spline, to_spline] = oneDquad_spline(x,knots);
    [from_spline, to_spline] = extend_to_nd(from_spline,to_spline,n_prior,n_post);
    from_spline(end+1,end+1) = 1;
    to_spline(end+1,end+1) = 1;
    n_splined = size(from_spline,2);
    [from_spline, to_spline] = projection_for_subset(from_spline,to_spline,0,n_g_red);

    % Reduce the decision vector
    [g1,psi,~,constant,g0] = change_basis(to_spline,from_spline,g1,psi,pi,constant,g0);
    pi = to_spline * pi * from_spline(1:n_v,1:n_splined);

else
    % Create identity matrix for code reuse below
    from_spline = speye(n_g_red + n_v);
    to_spline = speye(n_g_red + n_v);
    n_splined = n_v;
end


%% Step 4: Solve Linear Systems
t0 = tic;
fprintf('Solving linear system...\n');

% Note that I (SeHyoun) will probably swap out some parts of schur_solver in the
%    near future, so it might have a different interface. Since the codebase is
%    new, there will be interative updates to simplify interface and syntax.
%    Underlying math will stay the same, but interfact may change with updates,
%    so one should note the release number of the codebase one is using.

[G1, ~, impact, eu, F] = schur_solver(g0,g1,c,psi,pi,1,1,1,n_splined);
% fprintf('...Done!\n')
fprintf('Existence and uniqueness? %2.0f and %2.0f\n',eu);
% fprintf('Time to solve linear system: %2.4f seconds\n\n\n',toc(t0));


%% Step 5: Internal Consistency Check
% A different internal consistency check will be implemented and updated in the
%    future. This function should only be taken as a sanity check.

% impact1=impact(:,1);

if check_consistency
    g1 = -derivs(:,1:nVars);
    psi = -derivs(:,2*nVars+nEErrors+1:2*nVars+nEErrors+n_shocks);
    from_red = inv_state_red * from_spline;
    to_red = to_spline * state_red;
%    internal_consistency_check(G1,impact,n_g_red,from_red,to_red,g1,psi,F,n_v,n_g,T,steadystate,plotting,IRF,dt)
    [epsilon] = internal_consistency_check(G1,impact(:,1),n_g_red,from_red,to_red,g1,psi(:,1),F,n_v,n_g,1000,vars_SS,1,1, 0.07);

 end

%%