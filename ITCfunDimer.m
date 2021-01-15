function dQ = ITCfunDimer(X,p)

%% 
% functions to fit ITC data as below:
%		PP + 2L <--> 2P + 2L <--> 2PL
%
% Written by Olof Stenström 2020-06-11

% Order or variable in p
% p = [K1 K2 dH1 dH2 Qoff V0 Cs Cc Ex]


Vi = [p(9); X];

K1 = p(1)*1e-3;	%Binding Constant
K2 = p(2)*1e-3;	%Dimerisation constant
dH1 = p(3)*1e3;	%Enthalpy of binding
dH2 = p(4)*1e3;	%Enthalpy of dimerisation
Qoff = p(5)*1e-6;	%Offset from Zero
V0 = p(6);	%Cell volume
Cs = p(7);	%Syringe Concentration
Cc = p(8);	%Cell concentration



% Calculate initial concentrations
Ptot = Cc;
Ltot = 0;

P = -1/(4*K2)+sqrt((1/(4*K2))^2+Ptot/(2*K2));
PP = K2*P^2;

L = 0;
PL = 0;


% Calculate concentrations for all injections
for j = 1:length(Vi)
    Ptot(j+1) = Ptot(j)*V0/(V0+Vi(j));
    Ltot(j+1) = (Cs*Vi(j)+V0*Ltot(j))/(V0+Vi(j));
    syms x
    p = vpasolve(2*x^2*K2+x+x*Ltot(j+1)/(K1+x) == Ptot(j+1),x,Ptot(j+1));
    P(j+1) = double(p(3));
    L(j+1) = K1*Ltot(j+1)/(K1+P(j+1));
    PP(j+1) = P(j+1)^2*K2;
    PL(j+1) = P(j+1)*L(j+1)/K1;


end

% Calculate heat for each injection   
Q0 = V0*dH2*PP(1); % Start med att beräkna energin vid tiden 0. Den kommer ju inte vara 0 i detta fallet
dQ = [];
for i = 1:length(Vi)
    Q = V0*dH1*PL(i+1) + V0*dH2*PP(i+1);
    dQ(i) = Q-Q0 +(Vi(i)/V0)*(Q-Q0)/2 + Qoff;
    Q0 = Q;
end
dQ = dQ(2:end);