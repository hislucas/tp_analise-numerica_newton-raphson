% Lucas Ferreira Simões 2026020102
% 28 de outubro de 2026

% Trabalho prático de Análise Numérica
% Arquivo main.m
% Feito em Octave, compatível com Matlab, atendendo ao requisito de linguagem da proposta.

% A ideia é implementar Newton-Raphson para várias variáveis! Como fazer? A fórmula tá no livro e nas anotações.
% O Xande falou alguma coisa de Série de Mac-Laurim (não entendi a utilidade? Acho que é pra criar a matriz hessidiana -> derivadas de segunda ordem.
% Matriz Jacobiana para sistemas de equações não-lineares. Derivadas de primeira ordem.

% Precisei descobrir como fazer sistemas em Octave. E aparentemente tem diferenças de Matlab e Octave em álgebra linear. :. cogito passar para python.

%{
cp = 1.03 %[K J/(KgK)]
p = 331 + 45.6*w - 4.03*w^2; %[Pa]

%atmosférica para a entrada do compressor
p1 = 101000; %[Pa]
t1 = 25; %[°C]
Ec = w*cp*(t2-25)

%saida do compressor / entrada do queimador
p2 = p; %[Pa]
t2 = ?; %[°C]
% 8000 = w*cp*(t3-t2); []
% Como colocar essa eq no sistema? Depende de w. Sabendo cp (vide enunciado) e t3 e t2, dá pra fazer. 
% Calor específico = 1.03*10^3 [J/(KgK) -> Lê-se Quilograma Kelvin]
% Talvez seja assim que vamos descobrir t3. Sim, ficou t3 = 8000/(w*cp) + t2
%%% Na nova estrutura, os valores não vão ser calculados sequencialmente. Portanto, melhor deixar na forma mais fácil de derivar (original).

%saida do queimador / entrada da turbina:
w = 8.5019 + 0.02332*p + 0.48*10^(-4)*p^2 - 0.02644*t + 0.1849*10^(-4)*t^2 + 0.000121*p*t - 0.2736*10^(-6)*p^2*t - 0.1137*10^(-6)*p*t^2 + 0.2124*10^(-9)*p^2*t^2; %[kg/s]
p3 = p2; %[Pa]
t3 = ?; %[°C]

Et = 1727.5 - 10.06*p + 0.033033*p^2 - 7.4709*t + 0.003919*t^2 + 0.050921*p*t - 0.8525*10^(-4)*p^2*t - 0.2356*10^(-4)*p*t^2 + 0.4473*10^(-7)*p^2*t^2; %[kW]
Es = Et-Ec; %[W]

%saida da turbina para a atmosfera
t4 = t1;
p4 = p1; %[Pa]
%}
% Agora essas variáveis vão virar matriz.


X(:, 1) = [
	10; % X(1,1) = w inicial
	300; % X(2,1) = p inicial
	1000; % X(3,1) = Ec inicial
	250; % X(4,1) = t2 inicial
	1500; % X(5,1) = Es inicial
	800; % X(6,1) = t3 inicial
	2400; % X(7,1) = Et inicial
]


function F = calcular_F(X)
	F = zeros(7, 1); % reserva o lugar da coluna

	% p = 331 + 45.6w - 4.03w^2     % vale pra p2 e p3, que são iguais
	F(1) = X(2) - 331 - 45.6*X(1) + 4.03*X(1)^2;

	% Ec = 1020 - 0.383p + 0.00513p^2
	F(2) = X(3) - 1020 + 0.383*X(2) - 0.00513*X(2)^2;

	% w = f(p, t). Passando tudo para o lado de w (X1) e invertendo os sinais.
	F(3) = X(1) - 8.5019 - 0.02332*X(2) - 0.48e-4*X(2)^2 + 0.02644*X(6) - 0.1849e-4*X(6)^2 - 0.000121*X(2)*X(6) + 0.2736e-6*X(2)^2*X(6) + 0.1137e-6*X(2)*X(6)^2 - 0.2124e-9*X(2)^2*X(6)^2;

	% Et = f(p, t). Passando tudo para o lado de Et (X7) e invertendo os sinais.
	F(4) = X(7) - 1727.5 + 10.06*X(2) - 0.033033*X(2)^2 + 7.4709*X(6) - 0.003919*X(6)^2 - 0.050921*X(2)*X(6) + 0.8525e-4*X(2)^2*X(6) + 0.2356e-4*X(2)*X(6)^2 - 0.4473e-7*X(2)^2*X(6)^2;

	% Ec = w * cp * (t2-25), onde cp = 1.03
	F(5) = X(3) - X(1)*1.03*( X(4)-25 ); % F(5) = X(3) - X(1)*1.03*X(4) - X(1)*1.03*(-25) = X(3) - X(1)*1.03*X(4) + X(1)*1.03*(25)

	% 8000 = w * cp * (t-t2)
	F(6) = 8000 - X(1) * 1.03 * ( X(6) - X(4) ); % F(6) = 8000 - X(1)*1.03*X(6) + X(1)*1.03*X(4)

	% Et = Ec + Es
	F(7) = X(7) - X(3) - X(5);

end


% matriz jacobiana
function J = calcular_Jacobiana(X)
	J = zeros(7, 7);

	%1a equacao
	J(1, 1) = 0 - 0 - 45.6 + 8.06*X(1); % em função de x1
	J(1, 2) = 1 - 0 - 0 + 0; % em função de x2
	% o resto é 0, que sao as derivações em funcao das outras variáveis que não estão nessa equação

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




k = 1;
tol = 10^(-10);
erro = 1;
while (erro > tol)

	X_atual = X(:, k);

	F = calcular_F(X_atual);

	J = calcular_Jacobiana(X_atual);

	% aplicar lu às variáveis L, U e P (returns do lu)
	[L, U, P] = lu(J); %não é lupo igual à cueca
	%como aplicar? pesquisar


end





