function q1_dijkstra
% Q1_DIJKSTRA - Validacao da 1a Questao (Dijkstra, origem 0). MATLAB R2013b.
% Algoritmo: dijkstra() do material "algoritmos em MATLAB para Grafos".
% Conferencia extra: graphshortestpath (Bioinformatics Toolbox), se houver.
close all; clc;
prepara_pastas();
arqOut = '../out/q1_dijkstra.txt';
if exist(arqOut, 'file'), delete(arqOut); end
diary(arqOut);

% Grafo (vertice k da prova = indice k+1 no MATLAB)
nomes = {'0' '1' '2' '3' '4' '5'};
W = zeros(6);
W(1,3) = 20;  W(1,4) = 40;  W(1,5) = 30;   % 0->2, 0->3, 0->4
W(2,3) = 10;                               % 1->2
W(3,5) = 15;                               % 2->4
W(4,5) = 10;  W(4,6) = 80;                 % 3->4, 3->5
W(5,2) = 82;  W(5,6) = 20;                 % 4->1, 4->5
% Aresta 5->1 sem peso na legenda: nao incluida.

[dist, pred] = dijkstra_tab(W, 1, nomes);

fprintf('\n--- Resultado final (origem 0) ---\n');
fprintf('%-8s %-8s %-12s %s\n', 'Destino', 'Custo', 'Antecessor', 'Caminho');
for v = 1:6
    p = caminho(pred, 1, v);
    if pred(v) == 0, a = '--'; else a = nomes{pred(v)}; end
    fprintf('%-8s %-8g %-12s %s\n', nomes{v}, dist(v), a, strjoin(nomes(p), ' -> '));
end

esperado = [0 112 20 40 30 50];                    % resolucao manual
ok = isequal(dist, esperado);
if tem_bioinfo()
    dNat = graphshortestpath(sparse(W), 1, 'Method', 'Dijkstra');
    fprintf('\nDistancias (graphshortestpath): %s\n', mat2str(dNat));
    ok = ok && isequal(dNat, esperado);
else
    fprintf('\n(Bioinformatics Toolbox indisponivel: conferencia so com a resolucao manual)\n');
end
fprintf('Confere com a resolucao manual: %s\n', simnao(ok));

% Figura: arvore de caminhos minimos em destaque
D = false(6);
for v = 1:6
    if pred(v) ~= 0, D(pred(v), v) = true; end
end
X = [0 2 1  1   1 2];
Y = [0 0.7 1.2 -1.2 0 -0.7];
f = desenha_grafo(X, Y, nomes, W, true, D, [], ...
    'Q1 - Arvore de caminhos minimos (Dijkstra, origem 0)');
print(f, '-dpng', '-r200', '../figs/q1_dijkstra.png');
diary off;
end

% ===================== SUBFUNCOES =====================
function [dist, pred] = dijkstra_tab(W, s, nomes)
n = size(W, 1);
dist = inf(1, n);  pred = zeros(1, n);  visitado = false(1, n);
dist(s) = 0;
fprintf('Iter | u | '); fprintf('%-9s', nomes{:}); fprintf('\n');
for k = 1:n
    d = dist;  d(visitado) = inf;
    [m, u] = min(d);
    if isinf(m), break; end
    visitado(u) = true;
    for v = find(W(u, :))
        if ~visitado(v) && dist(u) + W(u, v) < dist(v)
            dist(v) = dist(u) + W(u, v);
            pred(v) = u;
        end
    end
    fprintf('%4d | %s | ', k, nomes{u});
    for v = 1:n
        if isinf(dist(v)), txt = 'inf';
        elseif pred(v) == 0, txt = sprintf('%g(-)', dist(v));
        else txt = sprintf('%g(%s)', dist(v), nomes{pred(v)}); end
        fprintf('%-9s', txt);
    end
    fprintf('\n');
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