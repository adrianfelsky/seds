function q2_bellman_ford
% Q2_BELLMAN_FORD - Validacao da 2a Questao (Bellman-Ford, origem A). MATLAB R2013b.
% Algoritmo: bellman_ford() do material + tabela por passadas.
% Conferencia extra: graphshortestpath 'Bellman-Ford' (Bioinformatics Toolbox).
close all; clc;
prepara_pastas();
arqOut = '../out/q2_bellman_ford.txt';
if exist(arqOut, 'file'), delete(arqOut); end
diary(arqOut);

nomes = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
L = 'ABCDEFG';
arestas = { 'A','B', 2;  'A','D', 5;  'A','F', 3;
            'B','E', 1;  'D','E', 1;  'F','B',-4;
            'E','C',-3;  'E','G', 3;
            'C','B', 7;  'C','G', 4;  'G','D',-1 };
W = zeros(7);
for e = 1:size(arestas, 1)
    W(L == arestas{e,1}, L == arestas{e,2}) = arestas{e,3};
end

% 1) Tabela por passadas (passada k usa os valores da passada k-1)
bellman_ford_tab(W, 1, nomes);

% 2) Implementacao do material
[dist, pred, cicloNeg] = bellman_ford(W, 1);
fprintf('\n--- Resultado final (bellman_ford do material) ---\n');
fprintf('%-6s %-6s %-10s %s\n', 'No', 'Custo', 'Anterior', 'Caminho');
for v = 1:7
    p = caminho(pred, 1, v);
    if pred(v) == 0, a = '--'; else a = nomes{pred(v)}; end
    fprintf('%-6s %-6g %-10s %s\n', nomes{v}, dist(v), a, strjoin(nomes(p), ' -> '));
end
fprintf('Ciclo negativo alcancavel: %s\n', simnao(cicloNeg));

esperado = [0 -1 -3 0 0 3 1];                      % resolucao manual
ok = isequal(dist, esperado);
if tem_bioinfo()
    dNat = graphshortestpath(sparse(W), 1, 'Method', 'Bellman-Ford');
    fprintf('\nDistancias (graphshortestpath Bellman-Ford): %s\n', mat2str(dNat));
    ok = ok && isequal(dNat, esperado);
else
    fprintf('\n(Bioinformatics Toolbox indisponivel: conferencia so com a resolucao manual)\n');
end
fprintf('Confere com a resolucao manual: %s\n', simnao(ok));

% Figura
D = false(7);
for v = 1:7
    if pred(v) ~= 0, D(pred(v), v) = true; end
end
%    A    B    C    D    E    F    G
X = [0    3    3    0    1.5  1.2  1.2];
Y = [2    2    0    0.6  1.2  3   -0.2];
f = desenha_grafo(X, Y, nomes, W, true, D, [], ...
    'Q2 - Arvore de caminhos minimos (Bellman-Ford, origem A)');
print(f, '-dpng', '-r200', '../figs/q2_bellman_ford.png');
diary off;
end

% ===================== SUBFUNCOES =====================
function bellman_ford_tab(W, s, nomes)
n = size(W, 1);
[I, J] = find(W);
peso = W(sub2ind(size(W), I, J));
dist = inf(1, n);  pred = zeros(1, n);  dist(s) = 0;
fprintf('Passada | '); fprintf('%-8s', nomes{:}); fprintf('\n');
imprime(0, dist, pred, nomes);
for k = 1:n-1
    dAnt = dist;
    for e = 1:numel(I)
        if dAnt(I(e)) + peso(e) < dist(J(e))
            dist(J(e)) = dAnt(I(e)) + peso(e);
            pred(J(e)) = I(e);
        end
    end
    imprime(k, dist, pred, nomes);
end
mudou = false;
for e = 1:numel(I)
    if dist(I(e)) + peso(e) < dist(J(e)), mudou = true; end
end
fprintf('Passada %d (verificacao): houve reducao? %s\n', n, simnao(mudou));
end

function imprime(k, dist, pred, nomes)
fprintf('%7d | ', k);
for v = 1:numel(dist)
    if isinf(dist(v)), txt = 'inf';
    elseif pred(v) == 0, txt = sprintf('%g(-)', dist(v));
    else txt = sprintf('%g(%s)', dist(v), nomes{pred(v)}); end
    fprintf('%-8s', txt);
end
fprintf('\n');
end

function [dist, pred, cicloNeg] = bellman_ford(W, s)
n = size(W, 1);
[I, J] = find(W);
peso = W(sub2ind(size(W), I, J));
dist = inf(1, n);  pred = zeros(1, n);  dist(s) = 0;
for k = 1:n-1
    mudou = false;
    for e = 1:numel(I)
        if dist(I(e)) + peso(e) < dist(J(e))
            dist(J(e)) = dist(I(e)) + peso(e);
            pred(J(e)) = I(e);
            mudou = true;
        end
    end
    if ~mudou, break; end
end
cicloNeg = false;
for e = 1:numel(I)
    if dist(I(e)) + peso(e) < dist(J(e))
        cicloNeg = true;  break;
    end
end
end

function p = caminho(pred, s, v)
p = v;
while v ~= s && pred(v) ~= 0
    v = pred(v);  p = [v p];
end
if p(1) ~= s, p = []; end
end

function ok = tem_bioinfo()
ok = license('test', 'Bioinformatics_Toolbox') && exist('graphshortestpath', 'file') > 0;
end

function s = simnao(x)
if x, s = 'SIM'; else s = 'NAO'; end
end

function prepara_pastas()
if ~exist('../figs', 'dir'), mkdir('../figs'); end
if ~exist('../out', 'dir'),  mkdir('../out');  end
end