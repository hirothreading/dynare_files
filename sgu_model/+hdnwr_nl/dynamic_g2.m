function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g2
%

if T_flag
    T = hdnwr_nl.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(60,1);
g2_j = zeros(60,1);
g2_v = zeros(60,1);

g2_i(1)=1;
g2_i(2)=1;
g2_i(3)=1;
g2_i(4)=2;
g2_i(5)=2;
g2_i(6)=2;
g2_i(7)=2;
g2_i(8)=2;
g2_i(9)=2;
g2_i(10)=2;
g2_i(11)=2;
g2_i(12)=2;
g2_i(13)=3;
g2_i(14)=3;
g2_i(15)=3;
g2_i(16)=4;
g2_i(17)=4;
g2_i(18)=4;
g2_i(19)=4;
g2_i(20)=4;
g2_i(21)=4;
g2_i(22)=4;
g2_i(23)=4;
g2_i(24)=5;
g2_i(25)=5;
g2_i(26)=5;
g2_i(27)=5;
g2_i(28)=5;
g2_i(29)=5;
g2_i(30)=5;
g2_i(31)=6;
g2_i(32)=6;
g2_i(33)=6;
g2_i(34)=6;
g2_i(35)=6;
g2_i(36)=6;
g2_i(37)=6;
g2_i(38)=6;
g2_i(39)=6;
g2_i(40)=6;
g2_i(41)=6;
g2_i(42)=6;
g2_i(43)=6;
g2_i(44)=6;
g2_i(45)=6;
g2_i(46)=6;
g2_i(47)=6;
g2_i(48)=6;
g2_i(49)=6;
g2_i(50)=6;
g2_i(51)=6;
g2_i(52)=6;
g2_i(53)=6;
g2_i(54)=7;
g2_i(55)=7;
g2_i(56)=8;
g2_i(57)=9;
g2_i(58)=9;
g2_i(59)=10;
g2_i(60)=10;
g2_j(1)=91;
g2_j(2)=98;
g2_j(3)=210;
g2_j(4)=73;
g2_j(5)=235;
g2_j(6)=230;
g2_j(7)=150;
g2_j(8)=236;
g2_j(9)=252;
g2_j(10)=151;
g2_j(11)=247;
g2_j(12)=253;
g2_j(13)=91;
g2_j(14)=98;
g2_j(15)=210;
g2_j(16)=73;
g2_j(17)=78;
g2_j(18)=158;
g2_j(19)=80;
g2_j(20)=192;
g2_j(21)=163;
g2_j(22)=165;
g2_j(23)=197;
g2_j(24)=1;
g2_j(25)=8;
g2_j(26)=120;
g2_j(27)=10;
g2_j(28)=154;
g2_j(29)=129;
g2_j(30)=161;
g2_j(31)=55;
g2_j(32)=56;
g2_j(33)=72;
g2_j(34)=57;
g2_j(35)=89;
g2_j(36)=52;
g2_j(37)=4;
g2_j(38)=61;
g2_j(39)=157;
g2_j(40)=62;
g2_j(41)=174;
g2_j(42)=73;
g2_j(43)=74;
g2_j(44)=90;
g2_j(45)=79;
g2_j(46)=175;
g2_j(47)=91;
g2_j(48)=96;
g2_j(49)=176;
g2_j(50)=10;
g2_j(51)=154;
g2_j(52)=163;
g2_j(53)=181;
g2_j(54)=55;
g2_j(55)=181;
g2_j(56)=55;
g2_j(57)=19;
g2_j(58)=199;
g2_j(59)=37;
g2_j(60)=217;
g2_v(1)=(-(y(13)*getPowerDeriv(y(6),params(4),2)));
g2_v(2)=(-T(38));
g2_v(3)=g2_v(2);
g2_v(4)=getPowerDeriv(y(5),(-params(2)),2);
g2_v(5)=(-(params(1)*(1+y(9))*getPowerDeriv(y(14),(-params(2)),2)/(1+y(15))));
g2_v(6)=(-(params(1)*T(37)/(1+y(15))));
g2_v(7)=g2_v(6);
g2_v(8)=(-((-(params(1)*(1+y(9))*T(37)))/((1+y(15))*(1+y(15)))));
g2_v(9)=g2_v(8);
g2_v(10)=(-((-(params(1)*T(6)))/((1+y(15))*(1+y(15)))));
g2_v(11)=g2_v(10);
g2_v(12)=(-((-((-(params(1)*(1+y(9))*T(6)))*(1+y(15)+1+y(15))))/((1+y(15))*(1+y(15))*(1+y(15))*(1+y(15)))));
g2_v(13)=y(13)*params(4)*getPowerDeriv(y(6),params(4)-1,2);
g2_v(14)=params(4)*T(39);
g2_v(15)=g2_v(14);
g2_v(16)=(-(y(12)*T(8)*1/(steady_state(2))*1/(steady_state(2))*getPowerDeriv(y(5)/(steady_state(2)),params(12),2)));
g2_v(17)=(-(y(12)*T(35)*T(41)));
g2_v(18)=g2_v(17);
g2_v(19)=(-(T(8)*T(35)));
g2_v(20)=g2_v(19);
g2_v(21)=(-(y(12)*T(9)*(1+params(9))/params(1)*1/(1+params(9))*1/(1+params(9))*getPowerDeriv((1+y(10))/(1+params(9)),params(11),2)));
g2_v(22)=(-(T(9)*T(41)));
g2_v(23)=g2_v(22);
g2_v(24)=(-((1+y(10))*(-((-y(8))*(y(1)+y(1))))/(y(1)*y(1)*y(1)*y(1))));
g2_v(25)=(-((1+y(10))*(-1)/(y(1)*y(1))));
g2_v(26)=g2_v(25);
g2_v(27)=(-((-y(8))/(y(1)*y(1))));
g2_v(28)=g2_v(27);
g2_v(29)=(-(1/y(1)));
g2_v(30)=g2_v(29);
g2_v(31)=T(11)*(T(25)*T(10)*T(22)*T(22)*T(44)+T(24)*T(24)*T(45));
g2_v(32)=T(24)*T(25)*T(36);
g2_v(33)=g2_v(32);
g2_v(34)=T(11)*(T(25)*T(22)*T(23)*1/(1-params(10))+T(24)*T(40)*T(45));
g2_v(35)=g2_v(34);
g2_v(36)=(-(T(1)*params(7)/(1+y(10))));
g2_v(37)=g2_v(36);
g2_v(38)=(-((-(y(1)*T(1)*params(7)))/((1+y(10))*(1+y(10)))));
g2_v(39)=g2_v(38);
g2_v(40)=T(11)*(T(25)*T(10)*(T(23)*(-(T(1)*params(7)))/((1+y(11))*(1+y(11)))+T(22)*T(42)*T(44))+T(24)*T(43)*T(45));
g2_v(41)=g2_v(40);
g2_v(42)=T(14)*getPowerDeriv(y(5),params(2),2);
g2_v(43)=T(36)*T(25)*T(40);
g2_v(44)=g2_v(43);
g2_v(45)=T(36)*T(25)*T(43);
g2_v(46)=g2_v(45);
g2_v(47)=T(11)*T(40)*T(40)*T(45);
g2_v(48)=T(11)*(T(40)*T(43)*T(45)+T(25)*1/(1-params(10))*T(23)*T(42));
g2_v(49)=g2_v(48);
g2_v(50)=(-((-T(2))/((1+y(10))*(1+y(10)))));
g2_v(51)=g2_v(50);
g2_v(52)=(-((-((-(T(2)*y(1)))*(1+y(10)+1+y(10))))/((1+y(10))*(1+y(10))*(1+y(10))*(1+y(10)))));
g2_v(53)=T(11)*(T(43)*T(43)*T(45)+T(25)*T(10)*(T(42)*T(42)*T(44)+T(23)*(-((-T(2))*(1+y(11)+1+y(11))))/((1+y(11))*(1+y(11))*(1+y(11))*(1+y(11)))));
g2_v(54)=(-(T(26)+T(26)+y(4)*T(1)*params(7)*T(1)*params(7)*getPowerDeriv(T(2),1-params(5),2)+T(3)*(-(params(7)*params(7)*getPowerDeriv(params(6)+params(7)*y(4),2-params(5),2)))));
g2_v(55)=getPowerDeriv(1+y(11),1-params(5),2);
g2_v(56)=(-(((y(4)+T(19))*(y(4)+T(19))*(T(34)*(1+T(31))+(y(4)+T(19))*(1-params(10))*(T(47)-(params(7)/(params(7)*(1-params(5)))*T(33)+params(7)/(params(7)*(1-params(5)))*T(33)+T(20)*(T(32)*T(46)+T(28)*T(28)*getPowerDeriv(T(16),1-params(5),2))))-(T(34)*(1+T(31))+(1-params(10))*(T(19)-T(20)*T(21))*T(47)))-((y(4)+T(19))*T(34)-(1-params(10))*(T(19)-T(20)*T(21))*(1+T(31)))*((y(4)+T(19))*(1+T(31))+(y(4)+T(19))*(1+T(31))))/((y(4)+T(19))*(y(4)+T(19))*(y(4)+T(19))*(y(4)+T(19)))));
g2_v(57)=(-(params(13)*(-1)/(y(2)*y(2))));
g2_v(58)=(-1)/(y(12)*y(12));
g2_v(59)=(-(params(14)*(-1)/(y(3)*y(3))));
g2_v(60)=(-1)/(y(13)*y(13));
g2 = sparse(g2_i,g2_j,g2_v,10,289);
end
