%Hansen's Simple RBC Model
%David Murakami
%%
var
y $y$ (long_name='output')
I $I$ (long_name='investment')
k $k$ (long_name='capital stock')
h $h$ (long_name='labour supply')
A $A$ (long_name='technology')
c $c$ (long_name='consumption')
r $r$ (long_name='real interest rate')
w $w$ (long_name='wage rate');

varexo eps $\epsilon$;

parameters alpha $\alpha$ (long_name='capital share')
beta $\beta$ (long_name='stochastic discount factor')
delta $\delta$ (long_name='depreciation')
rho $\rho$ (long_name='technology shock persistence')
eta $\eta$ (long_name='risk aversion coefficient')
a $a$ (long_name='labour disutility parameter')
sigmaeps $\sigma_{\epsilon}$ (long_name='volatility of shock');

alpha = 0.36;
beta =0.99;
delta = 0.025;
rho = 0.95;
eta = 1;
a = 2;
sigmaeps = 0.01;

%%
model;
1/c = beta*((1/c(+1))*(r(+1) +(1-delta))); %consumption euler equation
(1-alpha)*(y/h) = A/(1-h)*c; %labour first order condition
c = y +(1-delta)*k(-1) - k; %resource constraint
k = (1-delta)*k(-1) + I; %capital law of motion
y = A*k(-1)^(alpha)*h^(1-alpha); %production function
r = alpha*(y/k(-1)); %interest rate
w = (1-alpha)*(y/h); %wage rate
log(A) = rho*log(A(-1)) + eps; %shock process
end;

%%
initval;
A = 1;
h = (1+(a/(1-alpha))*(1-(beta*delta*alpha)/(1-beta*(1-delta))))^(-1);
k = h*((1/beta -(1-delta))/(alpha*A))^(1/(alpha-1));
I = delta*k;
y = A*k^(alpha)*h^(1-alpha);
c = y - delta*k;
r = 1/beta - 1 + delta;
w = (1-alpha)*(y/h);
end;

steady;

shocks;
var eps = sigmaeps^2;
end;

stoch_simul(order=1,irf=40);
write_latex_dynamic_model;
