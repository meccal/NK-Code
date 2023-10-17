%LUCA MECCA
%lmecca@london.edu
%Replicate the classical New Keynesian Model (Gali textbook 2015, chapter 3)
%PERTURBATION METHOD
%September 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define variables
var
    %variables needed to reproduce the figures
    y       ${y}$     (long_name='log output')
    infl    ${\pi}$   (long_name='log output')
    pi_ann  ${\pi^a}$ (long_name='log inflation (annualized)')
    n       ${n}$     (long_name='log employment')
    wp      ${n}$     (long_name='log real wage')
    i_ann   ${i}$     (long_name='log nominal interest rate (annualized)')
    r_ann   ${r}$     (long_name='log real interest rate (annualized)')


    Y   ${Y}$   (long_name='output')
    C   ${C}$   (long_name='consumption')
    PI  ${\Pi}$ (long_name='inflation')  
    PI_STAR  ${\Pi^*}$ (long_name='inflation for price-resetters')  
    N   ${L}$   (long_name='labour')
    I   ${I}$   (long_name='nominal interest rate')
    R   ${R}$   (long_name='real interest rate')
    MC_r  ${\Psi_t^r}$   (long_name='(real) marginal costs')
    WP  ${W}$   (long_name='real wage')
    X1  ${X_1}$ (long_name='First auxiliary variable')
    X2  ${X_2}$ (long_name='Second auxiliary variable')
    D   ${D_t}$ (long_name='Price dispersion')
    
    %exegenous processes
    a   ${a_t}$  (long_name='(log) Technology')
    z   ${z_t}$  (long_name='(log)Preference')
    v   ${v_t}$  (long_name='(log)Monetray Policy')
    
;

% Define exogenous variables
varexo
    eps_a    ${a_t}$   (long_name='Techonology shock')
    eps_z    ${z_t}$   (long_name='Preference shock')
    eps_v    ${v_t}$   (long_name='Monetrary policy shock')    
;


% Define parameters
parameters
    %Household
    SIGMA    ${\sigma}$    (long_name='CRRA coefficient')
    PHI      ${\phi}$      (long_name='Inverse of Frisch elasticity')
    BETA     ${\beta}$     (long_name='Subjective Discount Factor')
    EPSILON  ${\epsilon}$  (long_name='Elasticity of substitution between cosnumption goods')
    
    %Firm
    ALPHA    ${\alpha}$    (long_name='Degree of decrasing returns of labour')
    THETA    ${\theta}$    (long_name='Fraction of firms which cannot reoptimize')
    
    %Monetary Policy
    PHI_PI   ${\phi_{\pi}}$ (long_name='Reaction of monetary policy to inflation')
    PHI_Y    ${\phi_y}$     (long_name='Reaction of monetary policy to output')

    %Technology
    RHO_A    ${\rho^a}$    (long_name='AR(1) coefficient for technology')
    SIGMA_A  ${\sigma^a}$  (long_name='volatility for technology')
    %Preference
    RHO_Z    ${\rho^z}$    (long_name='AR(1) coefficient for preference shocks')
    SIGMA_Z  ${\sigma^z}$  (long_name='volatility for preference shocks')
    %Monetary Policy
    RHO_V    ${\rho^v}$    (long_name='AR(1) coefficient for monetary policy')
    SIGMA_V  ${\sigma^v}$  (long_name='volatility for monetary policy')
;

%Assign a value to the parameters
%Household
SIGMA=1; 
PHI=5;
BETA=0.99;
EPSILON=9;
ALPHA=0.25;
THETA=0.75;
PHI_PI=1.5;
PHI_Y=0.5/4;
RHO_A=0.9;
SIGMA_A=0.00025;
RHO_Z=0.5;
SIGMA_Z=0.00025;
RHO_V=0.5;
SIGMA_V=0.00025;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% MODEL %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;
    %Optimality conditions household
    [name='Euler equation']
    1=BETA*C^SIGMA*I/exp(z)*(exp(z(+1))/PI(+1)/C(+1)^SIGMA);
    [name='Intra-temporal consumption-leisure']
    N^PHI=WP/(C^SIGMA);
    
    %Optimality conditions firms
    [name='Real marginal cost']
    MC_r=WP/(exp(a)*N^(-ALPHA)*(1-ALPHA)*D^(ALPHA-1));
    [name='Law of motion first auxiliary variable']
    X1=exp(z)*Y*C^(-SIGMA)*MC_r  + BETA*THETA*(PI(+1)^(EPSILON/(1-ALPHA))*X1(+1));
    [name='Law of motion second auxiliary variable']
    X2=C^(-SIGMA)*exp(z)*Y + BETA*THETA*(PI(+1)^(EPSILON-1)*X2(+1));
    [name='Relationship between auxiliary variables']
    PI_STAR^((1-ALPHA+ALPHA*EPSILON)/(1-ALPHA))*X1=(EPSILON/(EPSILON-1))*X2;
    [name='Aggregate Price Dynamics']
    PI_STAR=((1-THETA*PI^(EPSILON-1))/(1-THETA))^(1/(1-EPSILON));
    [name='Law of motion of price dispersion']
    D=THETA*PI^(EPSILON/(1-ALPHA))*D(-1) + (1-THETA)*PI_STAR^(-EPSILON/(1-ALPHA));
    
    %Monetary Policy Rule
    [name='Interest Rate Rule']
    I=1/BETA*PI^(PHI_PI)*Y^(PHI_Y)*exp(v)/(steady_state(Y)^PHI_Y);
    [name='Real Interest Rate']
    R=I/PI(+1);
    
    %Market Clearing
    [name='Goods Market']
    Y=C;
    [name='Labour Market']
    N=(Y/exp(a))^(1/(1-ALPHA))*D;
    
    %Exogenous Processes
    [name='Stochastic Process for Technology']
    a=RHO_A*a(-1)+eps_a;
    [name='Stochastic Process for Preferences']
    z=RHO_Z*z(-1)+eps_z;
    [name='Stochastic Process for Monetary Policy']
    v=RHO_V*v(-1)+eps_v;
    
    %To rproduce the IRFs of Gali textbook
    [name='log output']
    y=log(Y);
    [name='log inflation']
    infl=log(PI);
    [name='log employment']
    pi_ann=log(PI)*4;
    [name='log employment']
    n=log(N);
    [name='log real wage']
    wp=log(WP);
    [name='log nominal interest rate (annualized)']
    i_ann=log(I)*4;
    [name='log real interest rate (annualized)']
    r_ann=log(R)*4;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% STEADY STATE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Steady State computation
%Initial values are taken from manual computation
 initval;
    Y=0.95;
    C=0.95;
    PI=1;
    PI_STAR=1; 
    N=0.93;  
    I=1/BETA;
    R=1/BETA;
    WP=0.68; 
    MC_r=0.68/(0.93^(-ALPHA)*(1-ALPHA));
    X1=3.45;  
    X2=3.88;  
    D=1;
    a=0;   
    z=0;  
    v=0;
    y=log(0.95);
    infl=0;
    pi_ann=0;
    n=log(0.93);
    wp=log(0.68);
    i_ann=0;
    r_ann=0;
  end;
  steady;

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% SIMULATIONS %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1. IRF to a monetary policy shock
%Define the volatility of the monetary policy shock
%size of the shock (epsilon) is 0.25
shocks;
var eps_v = SIGMA_V^2*0.25^2;
end;

set(groot, 'defaultLineLineWidth', 2.0)
stoch_simul(order=1, irf=16, periods=0, irf_plot_threshold=0) y pi_ann;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Policy Functions
shocks;
var eps_a = SIGMA_A^2;
var eps_z = SIGMA_Z^2;
var eps_v = SIGMA_V^2;
end;
set(groot, 'defaultLineLineWidth', 2.0)
stoch_simul(order=1, irf=0, periods=1000, nocorr, nodecomposition, nomoments, nograph) y infl;


%Note: columns in the ghx matrix are the state variables following in the
%DR-order. Because of this order, the variable z is last because it is the
%only mixed-varialbe (has one forward looking coefficient), while the other
%variables are purely backward variables (appear only at t and t+1)




