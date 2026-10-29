% Lucas Ferreira Simões 2026020102
% 28 de outubro de 2026

% Trabalho prático de Análise Numérica
% Arquivo main.m
% Feito em Octave, compatível com Matlab, atendendo ao requisito de linguagem da proposta.

% A ideia é implementar Newton-Raphson para várias variáveis! Como fazer? A fórmula tá no livro e nas anotações.
% O Xande falou alguma coisa de Série de Mac-Laurim (não entendi a utilidade? Acho que é pra criar a matriz hessidiana -> derivadas de segunda ordem.
% Matriz Jacobiana para sistemas de equações não-lineares. Derivadas de primeira ordem.


cp = 1.03*10^3 %[J/(KgK)]
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

%saida do queimador / entrada da turbina:
w = 8.5019 + 0.02332*p + 0.48*10^(-4)*p^2 - 0.02644*t + 0.1849*10^(-4)*t^2 + 0.000121*p*t - 0.2736*10^(-6)*p^2*t - 0.1137*10^(-6)*p*t^2 + 0.2124*10^(-9)*p^2*t^2; %[kg/s]
p3 = p2; %[Pa]
t3 = ?; %[°C]

Et = 1727.5 - 10.06*p + 0.033033*p^2 - 7.4709*t + 0.003919*t^2 + 0.050921*p*t - 0.8525*10^(-4)*p^2*t - 0.2356*10^(-4)*p*t^2 + 0.4473*10^(-7)*p^2*t^2; %[kW]
Es = Et-Ec; %[W]

%saida da turbina para a atmosfera
t4 = t1;
p4 = 101000; %[Pa]

