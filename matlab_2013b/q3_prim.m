function q3_prim
% Q3_PRIM - Validacao da 3a Questao (Prim, inicio em A). MATLAB R2013b.
% Algoritmo: prim() do material.
% Conferencia extra: graphminspantree 'Prim' (Bioinformatics Toolbox).
close all; clc;
prepara_pastas();
arqOut = '../out/q3_prim.txt';
if exist(arqOut, 'file'), delete(arqOut); end
diary(arqOut);

nomes = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
L = 'ABCDEFG';
arestas = { 'A','D', 5;  'A','B', 7;  'B','C', 8;  'B','D', 9;
            'B','E', 7;  'C','E', 5;  'D','E',15;  'D','F', 6;
            'E','F', 8;  'E','G', 9;  'F','G',11 };
P = zeros(7);
for e = 1:size(arestas, 1)
    P(L == arestas{e,1}, L == arestas{e,2}) = arestas{e,3};
end
P = P + P';

[ag, pesoTotal] = prim_tab(P, 1, nomes);
fprintf('\n--- AGM (ordem de inclusao) ---\n');
for k = 1:size(ag, 1)
    fprintf('%da: (%s, %s) = %g\n', k, nomes{ag(k,1)}, nomes{ag(k,2)}, P(ag(k,1), ag(k,2)));
end
fprintf('Peso total: %g\n', pesoTotal);

ok = (pesoTotal == 39);
if tem_bioinfo()
    T = graphminspantree(sparse(tril(P)), 1, 'Method', 'Prim');
    T = max(T, T');
    pNat = full(sum(nonzeros(tril(T))));
    fprintf('\nPeso total (graphminspantree Prim): %g\n', pNat);
    ok = ok && (pNat == 39);
else
    fprintf('\n(Bioinformatics Toolbox indisponivel: conferencia so com a resolucao manual)\n');
end
fprintf('Confere com a resolucao manual (39): %s\n', simnao(ok));

% Figura
D = false(7);
for k = 1:size(ag, 1), D(ag(k,1), ag(k,2)) = true; end
%    A    B    C    D    E    F    G
X = [0    1.5  3    0.5  2    1    2.8];
Y = [2    2.4  2    1    1    0    0  ];
f = desenha_grafo(X, Y, nomes, P, false, D, [], ...
    'Q3 - Arvore Geradora Minima (Prim, inicio A)');
print(f, '-dpng', '-r200', '../figs/q3_prim.png');
diary off;
end

% ===================== SUBFUNCOES =====================
function [arestas, pesoTotal] = prim_tab(W, inicio, nomes)
n = size(W, 1);
naArvore = false(1, n);  custo = inf(1, n);  pai = zeros(1, n);
custo(inicio) = 0;  arestas = zeros(0, 2);  pesoTotal = 0;
fprintf('Passo | u | '); fprintf('%-8s', nomes{:}); fprintf('\n');
for k = 1:n
    c = custo;  c(naArvore) = inf;
    [m, u] = min(c);
    if isinf(m), break; end
    naArvore(u) = true;
    if pai(u) ~= 0
        arestas(end+1, :) = [pai(u) u];
        pesoTotal = pesoTotal + W(pai(u), u);
    end
    for v = find(W(u, :))
        if ~naArvore(v) && W(u, v) < custo(v)
            custo(v) = W(u, v);
            pai(v) = u;
        end
    end
    fprintf('%5d | %s | ', k, nomes{u});
    for v = 1:n
        if naArvore(v) && v ~= u, txt = '--';
        elseif isinf(custo(v)), txt = 'inf';
        elseif pai(v) == 0, txt = sprintf('%g(-)', custo(v));
        else txt = sprintf('%g(%s)', custo(v), nomes{pai(v)}); end
        fprintf('%-8s', txt);
    end
    fprintf('\n');
end
end

function ok = tem_bioinfo()
ok = license('test', 'Bioinformatics_Toolbox') && exist('graphminspantree', 'file') > 0;
end

function s = simnao(x)
if x, s = 'SIM'; else s = 'NAO'; end
end

function prepara_pastas()
if ~exist('../figs', 'dir'), mkdir('../figs'); end
if ~exist('../out', 'dir'),  mkdir('../out');  end
end