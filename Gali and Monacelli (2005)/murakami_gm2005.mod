// Monetary Policy and the Small Open Economy
// Based on simple version of Gali and Monacelli (2005)
// Code written by David Murakami (Oxford, MPhil Economics)
// For use by Dynare 4.6.3
// Loops based on the replication by Johannes Pfeifer

@#define OPTIMAL = 0
@#define DITR = 0
@#define CITR = 0
@#define PEG = 1

% define string for saving different policies
@#if DITR == 1
    case_title='domestic inflation-based Taylor rule (DITR)';
@#else
    @#if CITR ==1
         case_title='CPI inflation-based Taylor rule (CITR)';
    @#else
        @#if PEG ==1
             case_title='exchange rate peg (PEG)';
        @#else
            @#if OPTIMAL ==1
                 case_title='Optimal Policy';
            @#else
                error('One case must be set to 1')
            @#endif
        @#endif
    @#endif
@#endif

// Define variables
var pih     $\pi^H$         (long_name='Domestic inflation')
    x       $x$             (long_name='Output gap')
    r       $i$             (long_name='Net nominal interest rate')
    pic     $\pi^C$         (long_name='CPI inflation')
    tot     $\tau$          (long_name='Terms of trade')
    ner     $e$             (long_name='Nominal exchange rate')
    der     $\Delta e$      (long_name='Nominal exchange rate depreciation')
    // Exogenous processes
    z       $z$             (long_name='Difference between foreign output and natural level of output')
    rn      $r^n$           (long_name='Net natural interest rate')
    rstar   $r^*$           (long_name='Net foreign interest rate')
    u       $u$             (long_name='Cost push shock')
    pistar  $\pi^*$         (long_name='Foreign inflation')
    ;

varexo epsz        $\varepsilon^{y^n}$     (long_name='Natural  rate of output shock')
       epsrn       $\varepsilon^{r^n}$     (long_name='Natural interest rate shock')
       epsrstar    $\varepsilon^{r^*}$     (long_name='Foreign interest rate shock')
       epsu        $\varepsilon^{u}$       (long_name='Cost push shock')
       ;

parameters BETA         $\beta$         (long_name='Discount factor')
           THETA        $\theta$        (long_name='Calvo parameter')
           VARPHI       $\varphi$       (long_name='Inverse Frisch elasticity of labour supply')
           ALPHA        $\alpha$        (long_name='Degree of openness')
           EPSILON      $\varepsilon$   (long_name='Elasticity of substition among varieties')
           PHIPI        $\phi_{\pi}$    (long_name='Taylor rule inflation coefficient')
           LAMBDA       $\lambda$       (long_name='Calvo factor (marginal cost coefficient in NKPC)')
           KAPPA        $\kappa$        (long_name='Slope of the NKPC')
           SIGMA        $\sigma$        (long_name='Risk aversion')
           OMEGA        $\Omega$        (long_name='Loss function adjustment parameter')
           LPI          $\lambda_{\pi}$ (long_name='Inflation weight in loss function')
           // AR coefficients
           RHOZ         $\rho_{z}$      (long_name='AR coefficient for gap between foreign and domestic output')
           RHORSTAR     $\rho_{r^*}$    (long_name='AR coefficient for foreign interest rate shock')
           RHORN        $\rho_{r^n}$    (long_name='AR coefficient for natural interest rate shock')
           RHOU         $\rho_{u}$      (long_name='AR coefficient for cost push shock')
           ;
// Set parameters
BETA = 0.99;
THETA = 0.75;
VARPHI = 3;
ALPHA = 0.4;
EPSILON = 6;
PHIPI = 1.5;
LAMBDA = (1-THETA)*(1-BETA*THETA)/THETA;
KAPPA = LAMBDA*(1+VARPHI);
SIGMA = 1;
OMEGA = (1-ALPHA)*(1+VARPHI);
LPI = EPSILON/(LAMBDA*(1+VARPHI));
RHOZ = 0.9;
RHORN = 0.9;
RHORSTAR = 0.9;
RHOU = 0.9;

model(linear);
[name='NKPC (for domestic inflation) (1)']
pih = KAPPA*x + BETA*pih(+1) + u;

[name='Dynamic IS equation (2)']
x = x(+1) - 1/SIGMA*(r - pih(+1) - rn);

[name='CPI inflation (3)']
pic = pih + ALPHA*(tot-tot(-1));

[name='Output gap (4)']
x = z + tot;

[name='pseudo UIP condition (5)']
der + pistar - pih + tot(-1) = rstar - pistar(+1) - r + pih(+1) + tot(+1);
/* Take a first difference of eq.(15) in GM's paper (i.e., derive tot-tot(-1)),
// then solve for tot:
tot - tot(-1) = der + pistar - pih
=> tot = der + pistar - pih + tot(-1)

Then use eq.(20) of GM and then set the two equations equal.
der + pistar - pih + tot(-1) = (rstar - pistar(+1)) - (r - pih(+1)) + tot(+1)
 */

[name='Nominal exchange rate depreciation (6)']
der = ner - ner(-1);

// Central bank rules
@#if DITR == 1
[name='Taylor rule for domestic inflation (7)']
r = PHIPI*pih; // domestic inflation-based Taylor rule (DITR)
@#else
    @#if CITR ==1
    [name='Taylor rule for CPI inflation (7)']
    r = PHIPI*pic; // CPI inflation-based Taylor rule (CITR)
    @#else
        @#if PEG ==1
        [name='Exchange rate rule (7)']
        ner=0; //exchange rate peg (PEG)
        @#else
            @#if OPTIMAL ==1
            [name='Optimal policy rule (7)']
            pih=0;
            @#else

            @#endif
        @#endif
    @#endif
@#endif

[name='Difference between foreign output and natural output (8)']
z = RHOZ*z(-1) + epsz;

[name='Natural rate of interest (9)']
rn = RHORN*rn(-1) + epsrn;

[name='Cost push shock (10)']
u = RHOU*u(-1) + epsu;

[name='Foreign interest rate (11)']
rstar = RHORSTAR*rstar(-1) + epsrstar;

[name='Foreign inflation rate (12)']
pistar = 0;

end;


steady;
check;
//model_diagnostics;


shocks;
var epsz = 1;
var epsrn = 1;
var epsu = 1;
var epsrstar = 1;
end;


write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;


stoch_simul(TeX,order=1,irf=20,irf_plot_threshold=0) pih x pic tot r der;


collect_latex_files;

var_string={'pih','x','pic','tot','r','der'};
