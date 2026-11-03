# Simulação de um sistema termodinâmico turbina–compressor–queimador

Trabalho Prático de Análise Numérica 2026/2 — Método de Newton–Raphson multivariável com fatoração LU.

| | |
|---|---|
| **Autor** | Lucas Ferreira Simões |
| **Matrícula** | 2026020102 |
| **Professor** | Alexandre Salles da Cunha |
| **Instituição** | Departamento de Ciência da Computação — Universidade Federal de Minas Gerais |
| **Linguagem** | Octave (compatível com MATLAB) |
| **Data** | Outubro de 2026 |

---

## Sobre o projeto

O sistema é composto por um **compressor**, um **queimador** e uma **turbina**. O ar à temperatura e pressão atmosféricas (`t1 = 25 °C`, `p1 = 101 kPa`) é captado pelo compressor e enviado ao queimador, onde são adicionados **8000 kW** de energia pela queima de combustível. Os gases de queima acionam a turbina, cujo eixo é acoplado ao do compressor: parte do trabalho gerado na turbina alimenta o compressor, e o restante (`Es`) fica disponível como potência livre no eixo.

O objetivo é determinar o **ponto de operação** do sistema, isto é, os valores de `w, p, Ec, t2, Es, t, Et` que satisfazem simultaneamente as sete equações do modelo (sistema não linear 7×7).

### Hipóteses simplificadoras

1. Ar atmosférico como gás perfeito, com calor específico constante `cp = 1,03 kJ/(kg·K)`.
2. Massa de combustível desconsiderada: a vazão `w` é constante ao longo do sistema.
3. Sem variação de pressão no queimador: `p2 = p3 = p`.
4. Sem perda de calor para o meio ambiente.

---

## Modelo matemático

**Variáveis:** `w` (vazão, kg/s), `p` (pressão de saída do compressor, kPa), `Ec` (potência consumida pelo compressor, kW), `t2` (temperatura na saída do compressor, °C), `Es` (potência livre no eixo, kW), `t = t3` (temperatura na entrada da turbina, °C) e `Et` (potência gerada na turbina, kW).

```text
(1)  p  = 331 + 45.6·w − 4.03·w²

(2)  Ec = 1020 − 0.383·p + 0.00513·p²

(3)  w  = 8.5019 + 0.02332·p + 0.48e-4·p² − 0.02644·t + 0.1849e-4·t²
          + 0.000121·p·t − 0.2736e-6·p²·t − 0.1137e-6·p·t² + 0.2124e-9·p²·t²

(4)  Et = 1727.5 − 10.06·p + 0.033033·p² − 7.4709·t + 0.003919·t²
          + 0.050921·p·t − 0.8525e-4·p²·t − 0.2356e-4·p·t² + 0.4473e-7·p²·t²

(5)  Ec = w·cp·(t2 − 25)          (balanço de energia no compressor)

(6)  8000 = w·cp·(t − t2)         (balanço de energia no queimador)

(7)  Et = Ec + Es                 (balanço de energia na turbina)
```

---

## Método numérico

O sistema é escrito na forma **F(x) = 0**, com `x = (w, p, Ec, t2, Es, t, Et)ᵀ`, e resolvido pelo **método de Newton–Raphson para várias variáveis**:

1. Calcula-se `F(x⁽ᵏ⁾)` e a matriz jacobiana `J(x⁽ᵏ⁾)`.
2. Resolve-se o sistema linear `J(x⁽ᵏ⁾) · Δ = −F(x⁽ᵏ⁾)` — **sem inverter a jacobiana** — usando fatoração **LU com pivoteamento parcial** (`P·J = L·U`).
3. Atualiza-se `x⁽ᵏ⁺¹⁾ = x⁽ᵏ⁾ + Δ`.
4. O processo para quando `‖Δ‖∞ ≤ 10⁻¹⁰`.

### Valores iniciais

| `w` | `p` | `Ec` | `t2` | `Es` | `t = t3` | `Et` |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 10 kg/s | 300 kPa | 1000 kW | 250 °C | 1500 kW | 800 °C | 2400 kW |

---

## Estrutura do repositório

```text
.
├── main.m                # Código-fonte do algoritmo (Newton–Raphson + LU)
├── main.tex              # Fonte LaTeX do relatório de resultados
├── resultados_latex.pdf  # Resultado da execução, com as iterações até a convergência
└── README.md
```

O arquivo `main.m` está organizado em:

- `calcular_F(x)` — avalia o vetor `F(x)` com as 7 equações do modelo;
- `calcular_Jacobiana(x)` — monta a matriz jacobiana 7×7 com as derivadas parciais analíticas;
- laço principal — iteração de Newton com `[L, U, P] = lu(J)` e `Δ = U \ (L \ (P * (−F)))`;
- impressão da tabela de iterações (`x`, erro e `Δ`) no terminal.

---

## Como executar

**Octave**

```bash
octave main.m
```

**MATLAB**

```matlab
>> main
```

A saída é uma tabela com `k`, o vetor `x⁽ᵏ⁾`, o erro `‖Δ‖∞` e o vetor de correção `Δ` de cada iteração.

---

## Resultados

Convergência em **5 iterações**, com erro final de `3,75 × 10⁻¹¹` (abaixo da tolerância de `1 × 10⁻¹⁰`).

### Iterações

Valores exibidos com 6 casas decimais; a tabela completa (10 casas) está em `resultados_latex.pdf`.

| k | w (kg/s) | p (kPa) | Ec (kW) | t2 (°C) | Es (kW) | t3 (°C) | Et (kW) | Erro `‖Δ‖∞` |
|:-:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 0 | 10.000000 | 300.000000 | 1000.000000 | 250.000000 | 1500.000000 | 800.000000 | 2400.000000 | — |
| 1 | 10.909359 | 352.172425 | 1507.404684 | 150.889385 | 1570.684982 | 877.573653 | 3078.089666 | 6,78e+02 |
| 2 | 10.767715 | 354.835615 | 1529.971225 | 162.793657 | 1597.745253 | 884.185351 | 3127.716477 | 4,96e+01 |
| 3 | 10.765296 | 354.854397 | 1530.068794 | 162.990006 | 1598.513701 | 884.474148 | 3128.582495 | 8,66e-01 |
| 4 | 10.765295 | 354.854401 | 1530.068809 | 162.990060 | 1598.514092 | 884.474268 | 3128.582900 | 4,05e-04 |
| 5 | 10.765295 | 354.854401 | 1530.068809 | 162.990060 | 1598.514092 | 884.474268 | 3128.582900 | 3,75e-11 |

<details>
<summary>Vetor de correção Δ por iteração</summary>

Ordem dos componentes: `(w, p, Ec, t2, Es, t3, Et)`. A iteração `k = 0` não possui vetor de correção.

| k | Δw | Δp | ΔEc | Δt2 | ΔEs | Δt3 | ΔEt |
|:-:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 9,09e-01 | 5,22e+01 | 5,07e+02 | -9,91e+01 | 7,07e+01 | 7,76e+01 | 6,78e+02 |
| 2 | -1,42e-01 | 2,66e+00 | 2,26e+01 | 1,19e+01 | 2,71e+01 | 6,61e+00 | 4,96e+01 |
| 3 | -2,42e-03 | 1,88e-02 | 9,76e-02 | 1,96e-01 | 7,68e-01 | 2,89e-01 | 8,66e-01 |
| 4 | -6,71e-07 | 4,06e-06 | 1,50e-05 | 5,41e-05 | 3,90e-04 | 1,20e-04 | 4,05e-04 |
| 5 | -2,82e-14 | -6,00e-13 | -1,84e-12 | 3,58e-12 | 3,75e-11 | 9,65e-12 | 3,59e-11 |

</details>

### Ponto de operação encontrado

| Variável | Valor (10 casas) | Valor arredondado |
|---|---:|---:|
| `w` | 10.7652949442 kg/s | 10,7653 kg/s |
| `p` | 354.8544012565 kPa | 354,8544 kPa |
| `Ec` | 1530.0688087663 kW | 1530,0688 kW |
| `t2` | 162.9900602383 °C | 162,9901 °C |
| `Es` | 1598.5140917183 kW | 1598,5141 kW |
| `t3` | 884.4742677339 °C | 884,4743 °C |
| `Et` | 3128.5829004846 kW | 3128,5829 kW |

Verificação rápida de consistência: `Ec + Es = 1530,0688 + 1598,5141 = 3128,5829 kW = Et`, conforme a Equação (7).

---

## Observações

- O erro cai de forma muito rápida nas últimas iterações (`4,96e+01 → 8,66e-01 → 4,05e-04 → 3,75e-11`), comportamento típico de convergência quadrática do método de Newton perto da solução.
- A variável `t3` foi a última a estabilizar, o que é esperado por ela aparecer em várias equações do sistema (3, 4 e 6).
- Em testes informais com outros valores iniciais para `t3`: partindo de 884, a variável voltou a 870 na iteração seguinte e depois convergiu normalmente; partindo de 870, foi para 871. Parece que o passo de subida de `t3` foi mais útil quando a variável começou longe da solução do que quando começou perto.
- O código foi desenvolvido e testado em **Octave** e depois executado no **MATLAB**, com resultados semelhantes. Os únicos avisos do MATLAB foram sobre vetores que mudam de tamanho dentro do laço (sem pré-alocação), sem impacto perceptível neste trabalho.

---

## Referências

- CAMPOS, Frederico Ferreira. *Algoritmos numéricos*. 4. ed. Edição do Autor, 2024. xiv, 428 p.
- [Matriz jacobiana — Wikipédia](https://pt.wikipedia.org/wiki/Matriz_jacobiana)
- [Método de Newton–Raphson — Wikipédia](https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Newton%E2%80%93Raphson)
- Notas das aulas de Análise Numérica.

README feito pelo Claude Sonnet 5.5 com pensamento extra (muito bom, inclusive).
