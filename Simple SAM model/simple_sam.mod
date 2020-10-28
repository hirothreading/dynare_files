%Simple Search and Match Model
%Author: David Murakami (University of Oxford)

// Declare variables
var y           $y$                 (long_name='output')
    w           $w$                 (long_name='wages')
    J           $J$                 (long_name='value of a job for the firm')
    V           $V$                 (long_name='job vacancy rate')
    U           $U$                 (long_name='unemployment')
    M           $M$                 (long_name='matching function (total hires)')
    thetta      $\theta$            (long_name='labour market tightness (V/U)')
    f           $f$                 (long_name='job finding probability')
    q           $q$                 (long_name='job filling probability')
;

// Declare shock process
varexo epsilon_y  $\epsilon_{y}$    (long_name='shock process for output')
;

// Declare parameters
parameters  betta    $\beta$         (long_name='discount factor')
            s        $s$             (long_name='separation rate')
            mu       $\mu$           (long_name='matching function elastcity')
            alppha   $\alpha$        (long_name='wage-output proportion')
            U_ss     $\bar{U}$       (long_name='steady state unemployment')
            q_ss     $\bar{q}$       (long_name='steady state job filling probability')
            y_ss     $\bar{y}$       (long_name='steady state output')
            rho_y    $\rho_{y}$      (long_name='output persistence')
            siggma_y $\sigma_{y}$    (long_name='standard deviation of output shock')
;

// Set parameter values
//load params;  // Can declare and load parameters from separate matlab file
//set_param_value(''     ,);
            betta    = 0.99;
            s        = 0.04;
            mu       = 0.5;
            alppha   = 2/3;
            U_ss     = 0.06;
            q_ss     = 0.7;
            y_ss     = 1;
            rho_y    = 0.95;
            siggma_y = 0.01;


model;
// Additional variables, parameters, and steady state values
// Turns out that a lot of this was redundant. Still, it's good to have anyway
// Steady state job finding probability
#f_ss = s/U_ss - s;

// Steady state labour market tightness
#thetta_ss = f_ss/q_ss;

// Steady state total hires
#M_ss = f_ss*U_ss;

// Steady state job vacancy rate
#V_ss = M_ss/q_ss;

// Specify value of matching function scalar parameter
#m = M_ss/(U_ss^mu*V_ss^(1-mu));

// Steady state wages
#w_ss = alppha*y_ss;

// Steady state job value for firm
#J_ss = (y_ss - w_ss)/(1-betta*(1-s));

// Specify cost of holding vacancy
#kapppa = q_ss*betta*J_ss;


// Specify observation equations
// Probability of finding a job (eq. 1)
f = M/U;

// Probability of filling a job (eq. 2)
q = M/V;

// Matching function (total hires) (eq. 3)
M = m*U^mu*V^(1-mu);

// Law of motion of unemployment (eq. 4)
U = (1-f)*U(-1) + s*(1-U(-1));

// Free entry condition (eq. 5)
kapppa = q*betta*J(+1);

// Worker output (eq. 6)
y = 1 - rho_y + rho_y*y(-1) + epsilon_y;

// Wages (eq. 7)
w = alppha*y;

// Employment transition equation (eq. 8)
J = y - w + betta*(1-s)*J(+1);

// Labour market tightness (eq. 9)
thetta = V/U;
end;


initval;
// We will have to manually input the numerical steady state values
y = 1.000000000;
U = 0.060000000;
q = 0.700000000;
f = 0.626666667;
thetta = 0.895238095;
M = 0.037600000;
V = 0.053714286;
w = 0.666666667;
J = 6.720430108;
end;


write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;


steady ;
check ;


shocks;
var epsilon_y ; stderr siggma_y ;
end;


options_.TeX = 1;
write_latex_dynamic_model;
write_latex_parameter_table;


stoch_simul(order=1,nomoments, irf=20, periods = 1000, nograph);


// Make plots
U_irf = U_ss - oo_.irfs.U_epsilon_y;
figure, plot(1:20,U_irf*100,'k','LineWidth',2),ylabel('percent')
