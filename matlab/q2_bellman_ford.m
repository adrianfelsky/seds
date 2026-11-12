% q2_bellman_ford.m - Validacao da 2a Questao (Bellman-Ford, origem A)
% Algoritmo: bellman_ford() do material "algoritmos em MATLAB para Grafos"
%            + versao por passadas (tabela k) + distances 'mixed' nativo.
% Requer MATLAB R2016b ou superior.
clear; close all; clc;
if ~exist('../figs','dir'), mkdir('../figs'); end
if ~exist('../out','dir'),  mkdir('../out');  end
arqOut = '../out/q2_bellman_ford.txt';
if exist(arqOut,'file'), delete(arqOut); end
diary(arqOut);

%% Grafo dirigido
nomes = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
L = 'ABCDEFG';
arestas = { 'A','B', 2;  'A','D', 5;  'A','F', 3;
            'B','E', 1;  'D','E', 1;  'F','B',-4;
            'E','C',-3;  'E','G', 3;
            'C','B', 7;  'C','G', 4;  'G','D',-1 };
W = zeros(7);
for e = 1:size(arestas,1)
    W(L == arestas{e,1}, L == arestas{e,2}) = arestas{e,3};
end

%% 1) Tabela por passadas (dist da passada k usa os valores da passada k-1)
bellman_ford_tab(W, 1, nomes);

%% 2) Implementacao do material
[dist, pred, cicloNeg] = bellman_ford(W, 1);
fprintf('\n--- Resultado final (bellman_ford do material) ---\n');
fprintf('%-6s %-6s %-10s %s\n','No','Custo','Anterior','Caminho');
for v = 1:7
    p = caminho(pred, 1, v);
    if pred(v) == 0, a = '--'; else, a = nomes{pred(v)}; end
    fprintf('%-6s %-6g %-10s %s\n', nomes{v}, dist(v), a, strjoin(nomes(p),' -> '));
end
fprintf('Ciclo negativo alcancavel: %s\n', simnao(cicloNeg));

%% 3) Funcao nativa (metodo 'mixed' aceita pesos negativos)
G = digraph(W, nomes);
dNat = distances(G, 'A', 'Method', 'mixed');
fprintf('\nDistancias (distances mixed nativo): %s\n', mat2str(dNat));
esperado = [0 -1 -3 0 0 3 1];                      % resolucao manual
fprintf('Confere com a resolucao manual: %s\n', simnao(isequal(dist, esperado) && isequal(dNat, esperado)));

%% Figura
figure('Name','Q2 - Bellman-Ford','Color','w');
Pl = plot(G,'EdgeLabel',G.Edges.Weight,'NodeFontSize',14,'NodeColor','k', ...
    'EdgeColor','b','LineWidth',1.5,'Marker','o','MarkerSize',8,'ArrowSize',12,'Layout','force');
TR = shortestpathtree(G, 'A', 'Method', 'mixed');
highlight(Pl, TR, 'EdgeColor','g', 'LineWidth', 4);
title('Q2 - Arvore de caminhos minimos (Bellman-Ford, origem A)');
salvar_fig(gcf, '../figs/q2_bellman_ford.png');
diary off;

%% ===================== FUNCOES =====================
function bellman_ford_tab(W, s, nomes)
    % Relaxamento por passadas: passada k = caminhos com ate k arestas
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
        else, txt = sprintf('%g(%s)', dist(v), nomes{pred(v)}); end
        fprintf('%-8s', txt);
    end
    fprintf('\n');
end

function [dist, pred, cicloNeg] = bellman_ford(W, s)
    % Menores caminhos com pesos possivelmente negativos. O(V*E).
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

function s = simnao(x)
    if x, s = 'SIM'; else, s = 'NAO'; end
end

function salvar_fig(f, arq)
    if exist('exportgraphics','file'), exportgraphics(f, arq, 'Resolution', 200);
    else, saveas(f, arq); end
end