% Autor: Lucas Ferreira Simões
% Número de matrícula: 2026020102
% 28 de outubro de 2026

% Trabalho prático de Análise Numérica
% Arquivo main.m
% Feito em Octave, compatível com Matlab, atendendo ao requisito de linguagem da proposta.

% A ideia é implementar Newton-Raphson para várias variáveis! Como fazer? A fórmula tá no livro e nas anotações.


X(:, 1) = [
	10; % X(1,1) = w inicial
	300; % X(2,1) = p inicial
	1000; % X(3,1) = Ec inicial
	250; % X(4,1) = t2 inicial
	1500; % X(5,1) = Es inicial
	800; % X(6,1) = t3 inicial -> referido por t
	2400; % X(7,1) = Et inicial
];


k = 1;
tol = 10^(-10); % tolerância máxima definida no enunciado
erro = 1;
while (erro > tol)

	X_atual = X(:, k);

	F = calcular_F(X_atual);

	J = calcular_Jacobiana(X_atual);

	% aplicar lu às variáveis L, U e P
	[L, U, P] = lu(J);

	% valor delta para ajuste
	delta = U \ ( L \ ( P * (-F) ) );

	% ajuste dos pesos por delta
	X(:, k+1) = X_atual + delta;

	%cálculo do erro
	erro = norm(delta, inf);

	% iteração
	k = k+1;
end


% display dos resultados. Depois vou colocar no latex e gerar o pdf






function F = calcular_F(X)
	F = zeros(7, 1); % reserva o lugar da coluna

	% p = 331 + 45.6w - 4.03w^2     % vale pra p2 e p3, que são iguais
	F(1) = X(2) - 331 - 45.6*X(1) + 4.03*X(1)^2;

	% Ec = 1020 - 0.383p + 0.00513p^2
	F(2) = X(3) - 1020 + 0.383*X(2) - 0.00513*X(2)^2;

	% w = 8.5019 + 0.02332*p + 0.48e-4*p^2 - 0.02644*t + 0.1849e-4*t^2 + 0.000121*p*t - 0.2736e-6*p^2*t - 0.1137e-6*p*t^2 + 0.2124e-9*p^2*t^2
	F(3) = X(1) - 8.5019 - 0.02332*X(2) - 0.48e-4*X(2)^2 + 0.02644*X(6) - 0.1849e-4*X(6)^2 - 0.000121*X(2)*X(6) + 0.2736e-6*X(2)^2*X(6) + 0.1137e-6*X(2)*X(6)^2 - 0.2124e-9*X(2)^2*X(6)^2;

	% Et = 1727.5 - 10.06*p + 0.033033*p^2 - 7.4709*t + 0.003919*t^2 + 0.050921*p*t - 0.8525e-4*p^2*t - 0.2356e-4*p*t^2 + 0.4473e-7*p^2*t^2
	F(4) = X(7) - 1727.5 + 10.06*X(2) - 0.033033*X(2)^2 + 7.4709*X(6) - 0.003919*X(6)^2 - 0.050921*X(2)*X(6) + 0.8525e-4*X(2)^2*X(6) + 0.2356e-4*X(2)*X(6)^2 - 0.4473e-7*X(2)^2*X(6)^2;

	% Ec = w * cp * (t2-25), onde cp = 1.03
	F(5) = X(3) - X(1)*1.03*( X(4)-25 ); % F(5) = X(3) - X(1)*1.03*X(4) - X(1)*1.03*(-25) = X(3) - X(1)*1.03*X(4) + X(1)*1.03*(25)

	% 8000 = w * cp * (t-t2)
	F(6) = 8000 - X(1) * 1.03 * ( X(6) - X(4) ); % F(6) = 8000 - X(1)*1.03*X(6) + X(1)*1.03*X(4)

	% Et = Ec + Es
	F(7) = X(7) - X(3) - X(5);

end


% matriz jacobiana :. derivar tudinho em função de cada variável
function J = calcular_Jacobiana(X)
	J = zeros(7, 7); % reserva a matriz

	%1a equacao
	J(1, 1) = 0 - 0 - 45.6 + 8.06*X(1); % em função de x1
	J(1, 2) = 1 - 0 - 0 + 0; % em função de x2
	% o resto é 0, que são as derivações em função das outras variáveis que não estão nessa equação

	%2a equacao
	J(2, 2) = 0 - 0 + 0.383 - 2*0.00513*X(2); %x2
	J(2, 3) = 1 - 0 + 0 - 0; %x3

	%3a equacao
	J(3, 1) = 1 - 0 - 0 - 0 + 0 - 0 - 0 + 0 + 0 - 0; %x1
	J(3, 2) = 0 - 0.02332 - 2*0.48e-4*X(2) + 0 - 0 - 0.000121*X(6) + 2*0.2736e-6*X(2)*X(6) + 0.1137e-6*X(6)^2 - 2*0.2124e-9*X(2)*X(6)^2; %x2
	J(3, 6) = 0 - 0 - 0 + 0.02644 - 2*0.1849e-4*X(6) - 0.000121*X(2) + 0.2736e-6*X(2)^2 + 2*0.1137e-6*X(2)*X(6) - 2*0.2124e-9*X(2)^2*X(6);  %x6

	%4a equacao
	J(4, 2) = 0 - 0 + 10.06 - 2*0.033033*X(2) + 0 - 0 - 0.050921*X(6) + 2*0.8525e-4*X(2)*X(6) + 0.2356e-4*X(6)^2 - 2*0.4473e-7*X(2)*X(6)^2; %x2
	J(4, 6) = 0 - 0 + 0 - 0 + 7.4709 - 2*0.003919*X(6) - 0.050921*X(2) + 0.8525e-4*X(2)^2 + 2*0.2356e-4*X(2)*X(6) - 2*0.4473e-7*X(2)^2*X(6); %x6
	J(4, 7) = 1 - 0 + 0 - 0 + 0 - 0 - 0 + 0 + 0 - 0; %x7

	%5a equacao
	J(5, 1) = 0 - 1.03*( X(4)-25 ); %x1
	J(5, 3) = 1 - 0; %x3
	J(5, 4) = 0 - X(1)*1.03 + 0; %x4

	%6a equacao
	J(6, 1) = 0 - 1.03*X(6) + 1.03*X(4); %x1
	J(6, 4) = 0 - 0 + X(1)*1.03; %x4
	J(6, 6) = 0 - X(1)*1.03; %x6

	%7a equacao
	J(7, 3) = 0 - 1 - 0; %x3
	J(7, 5) = 0 - 0 - 1; %x5
	J(7, 7) = 1 - 0 - 0; %x7

end
