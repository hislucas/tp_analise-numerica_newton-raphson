% Autor: Lucas Ferreira Simões
% Número de matrícula: 2026020102
% 28 de outubro de 2026

% Trabalho prático de Análise Numérica
% Arquivo main.m
% Feito em Octave, compatível com Matlab, atendendo ao requisito de linguagem da proposta.

% A ideia é implementar Newton-Raphson para várias variáveis! Como fazer? A fórmula tá no livro e nas anotações.

clear lucas_X lucas_X_original lucas_F lucas_J lucas_tol lucas_iteracoes lucas_k lucas_erro lucas_L lucas_U lucas_P lucas_delta

lucas_X(:, 1) = [
	10; % X(1,1) = w inicial
	300; % X(2,1) = p inicial
	1000; % X(3,1) = Ec inicial
	250; % X(4,1) = t2 inicial
	1500; % X(5,1) = Es inicial
	800; % X(6,1) = t3 inicial -> referido por t
	2400; % X(7,1) = Et inicial
];
lucas_X_original = lucas_X;


function lucas_F = calcular_F(lucas_X)
	lucas_F = zeros(7, 1); % reserva o lugar da coluna

	% p = 331 + 45.6w - 4.03w^2     % vale pra p2 e p3, que são iguais
	lucas_F(1) = lucas_X(2) - 331 - 45.6*lucas_X(1) + 4.03*lucas_X(1)^2;

	% Ec = 1020 - 0.383p + 0.00513p^2
	lucas_F(2) = lucas_X(3) - 1020 + 0.383*lucas_X(2) - 0.00513*lucas_X(2)^2;

	% w = 8.5019 + 0.02332*p + 0.48e-4*p^2 - 0.02644*t + 0.1849e-4*t^2 + 0.000121*p*t - 0.2736e-6*p^2*t - 0.1137e-6*p*t^2 + 0.2124e-9*p^2*t^2
	lucas_F(3) = lucas_X(1) - 8.5019 - 0.02332*lucas_X(2) - 0.48e-4*lucas_X(2)^2 + 0.02644*lucas_X(6) - 0.1849e-4*lucas_X(6)^2 - 0.000121*lucas_X(2)*lucas_X(6) + 0.2736e-6*lucas_X(2)^2*lucas_X(6) + 0.1137e-6*lucas_X(2)*lucas_X(6)^2 - 0.2124e-9*lucas_X(2)^2*lucas_X(6)^2;

	% Et = 1727.5 - 10.06*p + 0.033033*p^2 - 7.4709*t + 0.003919*t^2 + 0.050921*p*t - 0.8525e-4*p^2*t - 0.2356e-4*p*t^2 + 0.4473e-7*p^2*t^2
	lucas_F(4) = lucas_X(7) - 1727.5 + 10.06*lucas_X(2) - 0.033033*lucas_X(2)^2 + 7.4709*lucas_X(6) - 0.003919*lucas_X(6)^2 - 0.050921*lucas_X(2)*lucas_X(6) + 0.8525e-4*lucas_X(2)^2*lucas_X(6) + 0.2356e-4*lucas_X(2)*lucas_X(6)^2 - 0.4473e-7*lucas_X(2)^2*lucas_X(6)^2;

	% Ec = w * cp * (t2-25), onde cp = 1.03
	lucas_F(5) = lucas_X(3) - lucas_X(1)*1.03*( lucas_X(4)-25 ); % F(5) = X(3) - X(1)*1.03*X(4) - X(1)*1.03*(-25) = X(3) - X(1)*1.03*X(4) + X(1)*1.03*(25)

	% 8000 = w * cp * (t-t2)
	lucas_F(6) = 8000 - lucas_X(1) * 1.03 * ( lucas_X(6) - lucas_X(4) ); % F(6) = 8000 - X(1)*1.03*X(6) + X(1)*1.03*X(4)

	% Et = Ec + Es
	lucas_F(7) = lucas_X(7) - lucas_X(3) - lucas_X(5);
end


% matriz jacobiana :. derivar tudinho em função de cada variável
function lucas_J = calcular_Jacobiana(lucas_X)
	lucas_J = zeros(7, 7); % reserva a matriz

	%1a equacao
	lucas_J(1, 1) = 0 - 0 - 45.6 + 8.06*lucas_X(1); % em função de x1
	lucas_J(1, 2) = 1 - 0 - 0 + 0; % em função de x2
	% o resto é 0, que são as derivações em função das outras variáveis que não estão nessa equação

	%2a equacao
	lucas_J(2, 2) = 0 - 0 + 0.383 - 2*0.00513*lucas_X(2); %x2
	lucas_J(2, 3) = 1 - 0 + 0 - 0; %x3

	%3a equacao
	lucas_J(3, 1) = 1 - 0 - 0 - 0 + 0 - 0 - 0 + 0 + 0 - 0; %x1
	lucas_J(3, 2) = 0 - 0.02332 - 2*0.48e-4*lucas_X(2) + 0 - 0 - 0.000121*lucas_X(6) + 2*0.2736e-6*lucas_X(2)*lucas_X(6) + 0.1137e-6*lucas_X(6)^2 - 2*0.2124e-9*lucas_X(2)*lucas_X(6)^2; %x2
	lucas_J(3, 6) = 0 - 0 - 0 + 0.02644 - 2*0.1849e-4*lucas_X(6) - 0.000121*lucas_X(2) + 0.2736e-6*lucas_X(2)^2 + 2*0.1137e-6*lucas_X(2)*lucas_X(6) - 2*0.2124e-9*lucas_X(2)^2*lucas_X(6);  %x6

	%4a equacao
	lucas_J(4, 2) = 0 - 0 + 10.06 - 2*0.033033*lucas_X(2) + 0 - 0 - 0.050921*lucas_X(6) + 2*0.8525e-4*lucas_X(2)*lucas_X(6) + 0.2356e-4*lucas_X(6)^2 - 2*0.4473e-7*lucas_X(2)*lucas_X(6)^2; %x2
	lucas_J(4, 6) = 0 - 0 + 0 - 0 + 7.4709 - 2*0.003919*lucas_X(6) - 0.050921*lucas_X(2) + 0.8525e-4*lucas_X(2)^2 + 2*0.2356e-4*lucas_X(2)*lucas_X(6) - 2*0.4473e-7*lucas_X(2)^2*lucas_X(6); %x6
	lucas_J(4, 7) = 1 - 0 + 0 - 0 + 0 - 0 - 0 + 0 + 0 - 0; %x7

	%5a equacao
	lucas_J(5, 1) = 0 - 1.03*( lucas_X(4)-25 ); %x1
	lucas_J(5, 3) = 1 - 0; %x3
	lucas_J(5, 4) = 0 - lucas_X(1)*1.03 + 0; %x4

	%6a equacao
	lucas_J(6, 1) = 0 - 1.03*lucas_X(6) + 1.03*lucas_X(4); %x1
	lucas_J(6, 4) = 0 - 0 + lucas_X(1)*1.03; %x4
	lucas_J(6, 6) = 0 - lucas_X(1)*1.03; %x6

	%7a equacao
	lucas_J(7, 3) = 0 - 1 - 0; %x3
	lucas_J(7, 5) = 0 - 0 - 1; %x5
	lucas_J(7, 7) = 1 - 0 - 0; %x7

end


lucas_k = 1;
lucas_iteracoes(1) = 0;
lucas_tol = 10^(-10); % tolerância máxima definida no enunciado
lucas_erro(1) = 1;
while (lucas_erro > lucas_tol)

	lucas_X_atual = lucas_X(:, lucas_k);

	lucas_F = calcular_F(lucas_X_atual);

	lucas_J = calcular_Jacobiana(lucas_X_atual);

	% aplicar lu às variáveis L, U e P
	[lucas_L, lucas_U, lucas_P] = lu(lucas_J);

	% valor delta para ajuste
	lucas_delta(:, lucas_k) = lucas_U \ ( lucas_L \ ( lucas_P * (-lucas_F) ) );

	% ajuste dos pesos por delta
	lucas_X(:, lucas_k+1) = lucas_X_atual + lucas_delta(:, lucas_k);

	%cálculo do erro
	lucas_erro(lucas_k+1) = norm(lucas_delta(:, lucas_k), inf);

	% iteração
	lucas_k = lucas_k+1;
	lucas_iteracoes(lucas_k) = lucas_k-1;

end


% display dos resultados. Depois vou colocar no latex e gerar o pdf
clear lucas_k;

fprintf("k | w             | p              | Ec              | t2             | Es              | t3             | Et              |  Erro     | Delta\n");
for lucas_k = 1:size(lucas_iteracoes, 2)
    fprintf("%i ", lucas_k-1);
    fprintf("| %.10f ", lucas_X(:, lucas_k));
    fprintf("| %9.2e ", lucas_erro(lucas_k));

    if (lucas_k>1)
        fprintf("| ");
        fprintf("%10.2e ", lucas_delta(:, lucas_k-1));
    end

    fprintf("| \n");

end

% Referências:
% https://pt.wikipedia.org/wiki/Matriz_jacobiana
% https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Newton%E2%80%93Raphson
% CAMPOS, Frederico Ferreira. Algoritmos numéricos. 4.ed. Edição do Autor, 2024. xiv, 428p.
% Aulas de Análise Numérica
