function q5_kruskal
% Q5_KRUSKAL - Validacao da 5a Questao (Kruskal). MATLAB R2013b.
% Algoritmo: Kruskal com Union-Find conforme pseudocodigo do material (1.pdf).
% Conferencia extra: graphminspantree 'Kruskal' (Bioinformatics Toolbox).
% Rotulos: A (sup. esq.), B (sup. centro), C (sup. dir.),
%          D (inf. esq.), E (base), F (inf. dir.)
close all; clc;
prepara_pastas();
arqOut = '../out/q5_kruskal.txt';
if exist(arqOut, 'file'), delete(arqOut); end
diary(arqOut);

nomes = {'A' 'B' 'C' 'D' 'E' 'F'};
L = 'ABCDEF';
arestas = { 'A','B',2;  'B','C',3;  'A','D',5;  'A','E',1;  'B','E',7;
            'C','E',3;  'C','F',9;  'D','E',4;  'E','F',4 };
W = zeros(6);
for e = 1:size(arestas, 1)
    W(L == arestas{e,1}, L == arestas{e,2}) = arestas{e,3};
end
W = W + W';

[ag, pesoTotal] = kruskal_tab(W, nomes);
fprintf('\nArestas da AGM: ');
for k = 1:size(ag, 1)
    fprintf('(%s,%s) ', nomes{ag(k,1)}, nomes{ag(k,2)});
end
fprintf('\nPeso total: %g\n', pesoTotal);

ok = (pesoTotal == 14);
if tem_bioinfo()
    T = graphminspantree(sparse(tril(W)), 'Method', 'Kruskal');
    T = max(T, T');
    pNat = full(sum(nonzeros(tril(T))));
    fprintf('Peso total (graphminspantree Kruskal): %g\n', pNat);
    ok = ok && (pNat == 14);
else
    fprintf('(Bioinformatics Toolbox indisponivel: conferencia so com a resolucao manual)\n');
end
fprintf('Confere com a resolucao manual (14): %s\n', simnao(ok));

% Figura (posicoes semelhantes ao desenho da prova)
D = false(6);
for k = 1:size(ag, 1), D(ag(k,1), ag(k,2)) = true; end
X = [1.5 3.0 4.3 1.2 2.8 4.1];
Y = [3.0 3.2 3.0 1.5 1.0 1.5];
f = desenha_grafo(X, Y, nomes, W, false, D, [], ...
    'Q5 - Arvore Geradora Minima (Kruskal)');
print(f, '-dpng', '-r200', '../figs/q5_kruskal.png');
diary off;
end

% ===================== SUBFUNCOES =====================
function [A, pesoTotal] = kruskal_tab(W, nomes)
n = size(W, 1);
[u, v] = find(triu(W));
w = W(sub2ind(size(W), u, v));
E = sortrows([w u v], [1 2 3]);       % peso crescente (empate: ordem alfabetica)
pai = 1:n;                            % Criar-Conjunto(v)
A = zeros(0, 2);  pesoTotal = 0;
fprintf('Passo | Aresta | Peso | Decisao\n');
for k = 1:size(E, 1)
    a = E(k, 2);  b = E(k, 3);
    ra = buscar(pai, a);  rb = buscar(pai, b);
    if ra ~= rb
        A(end+1, :) = [a b]; %#ok<AGROW>
        pesoTotal = pesoTotal + E(k, 1);
        pai(rb) = ra;                 % Unir-Conjuntos(u, v)
        dec = 'Aceita';
    else
        dec = 'Rejeitada (ciclo)';
    end
    fprintf('%5d | (%s,%s)  | %4g | %s\n', k, nomes{a}, nomes{b}, E(k, 1), dec);
end
end

function r = buscar(pai, x)
% Buscar-Conjunto(x)
r = x;
while pai(r) ~= r, r = pai(r); end
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