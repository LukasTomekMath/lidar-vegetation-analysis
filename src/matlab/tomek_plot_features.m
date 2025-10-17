%% features plot
n1 = 22;
n2 = 11;
N = n1 + n2;

C1start = 1;
C2start = n1 + 1;

C1end = n1;
C2end = n1 + n2;

ind = [C1start, C1end; C2start, C2end];

X = dataAll(:,1);
Y = dataAll(:,2);
Z = dataAll(:,3);

figure
hold on

scatter(X(ind(1,1):ind(1,2)), Y(ind(1,1):ind(1,2)), 20, "red", "filled");
scatter(X(ind(2,1):ind(2,2)), Y(ind(2,1):ind(2,2)), 50, "black", "*");
hold off
% axis([-0.1 1.2 -0.1 1.2])
legend('91E0', 'Monokultura', 'Location','southeast')


%% 3D graf
K = coeff(1,2).const;
L = coeff(1,2).linear; 
f = @(x,y,z) K + L(1)*x + L(2)*y + L(3)*z;
figure
fimplicit3(f)
hold on
scatter3(X(ind(1,1):ind(1,2)), Y(ind(1,1):ind(1,2)), Z(ind(1,1):ind(1,2)), 20, "red", "filled");
scatter3(X(ind(2,1):ind(2,2)), Y(ind(2,1):ind(2,2)), Z(ind(2,1):ind(2,2)), 50, "black", "*");
hold off
legend('91E0', 'Monokultura', 'Location','southeast')

%% parallel plot

figure; 
hold on
% features = ;

features = 88;

plot([0,0],[0,0],'-r');
plot([0,0],[0,0],'-k');

plot(1:1:features, dataScaled(ind(1,1):ind(1,2),1:features),'xr')

plot(1:1:features, dataScaled(ind(2,1):ind(2,2),1:features),'xk')

legend('91E0', 'Monokultura', 'Location','southeast')
hold off