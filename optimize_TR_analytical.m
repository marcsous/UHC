%% code to optimize TR for T1int = T1null
clear

syms T1int T1null TElast TR a b c positive

% TI to ensure signal at T1null is zero
TI = T1null*(log(sym(2))-log(sym(1)+exp(-(TR-TElast)/T1null)));

% signal at T1int
S = sym(1) - sym(2)*exp(-TI/T1int) + exp(-(TR-TElast)/T1int);

% efficiency at T1int
E = S/sqrt(TR);

% taylor expansion about T1int = T1null
E = taylor(E,T1int,T1null,'order',2);

% isolate linear coefficient
f = E / (T1null-T1int);

% maximum requires df/dTR = 0
df = diff(f,TR);

% pull out common terms to help the computer
g = df * sym(2) * exp((TR-TElast)/T1null) * T1null * TR^(sym(3)/sym(2));

% result is a function of a = TR/T1null and c = TElast/T1null only
g = simplify(subs(g,[TR TElast],[a*T1null c*T1null]))

% exact solution at c=0
aopt = vpasolve(subs(g,c,0),a,4)

% for c->0 assume nearby solution: a = aopt + b*c
gapprox = taylor(subs(g,a,aopt+b*c),c,0,'order',2);

% discard numerically "zero" constant term
threshold = 10^(vpa(-digits/2)); % sqrt(eps)
gapprox = mapSymType(gapprox,'vpa',@(x)logical(abs(x)>threshold)*x)

% solve for b
bopt = solve(gapprox==0,b)

%% numerical validation
small_c = 0:0.1:1;

% exact and approximate solutions
for j = 1:numel(small_c)
    approx(j) = aopt + bopt*small_c(j);    
    exact(j) = vpasolve(subs(g,c,small_c(j)),a,approx(j));
end

% display
subplot(1,2,2);
plot(small_c,[approx;exact]);
xlabel('TE_{last} / T_1^{null}');
ylabel('TR / T_1^{null}');
grid on
xticks(0.25*(0:4));
title('Compare optimal TR estimates');
lg2 = legend({'Approximate','Exact'},'Location','NorthWest');
lg2.Title.String = 'Equation';

%% display efficiency curve
subplot(1,2,1)

TR = 1:5000; % ms
T1null = [505;721]; % ms
TElast = 0; % ms

func = matlabFunction(f);
fval = func(1e-3*T1null,1e-3*TElast,1e-3*TR); % s

plot(TR,fval);
grid on;
xlabel('TR (ms)');
ylabel('f (TR, T_1^{null}, TE_{last})');
yticks(0.25*(0:4));
xlim([0 max(TR)]);
ylim([0 1]);
title('Efficiency vs TR');
lg1 = legend(num2str(T1null,'%.0f'));
lg1.Title.String = 'T_1^{null} (ms)';
