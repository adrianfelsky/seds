% q4_ford_fulkerson.m - Validacao da 4a Questao (Ford-Fulkerson, s = 0, t = 5)
% Algoritmo: ford_fulkerson() do material "algoritmos em MATLAB para Grafos"
%            + verificacao dos caminhos da resolucao manual + maxflow nativo.
% Requer MATLAB R2016b ou superior.
clear; close all; clc;
if ~exist('../figs','dir'), mkdir('../figs'); end
if ~exist('../out','dir'),  mkdir('../out');  end
arqOut = '../out/q4_ford_fulkerson.txt';
if exist(arqOut,'file'), delete(arqOut); end
diary(arqOut);

%% Rede de fluxo (vertice k da prova = indice k+1 no MATLAB)
nomes = {'0' '1' '2' '3' '4' '5'};
C = zeros(6);
C(1,2) = 16;  C(1,3) = 13;     % 0->1, 0->2
C(2,3) = 10;  C(3,2) = 4;      % 1->2, 2->1
C(2,4) = 12;                   % 1->3
C(4,3) = 9;                    % 3->2
C(3,5) = 14;                   % 2->4
C(5,4) = 7;                    % 4->3
C(4,6) = 20;  C(5,6) = 4;      % 3->5, 4->5
s = 1;  t = 6;

%% 1) Caminhos aumentantes da resolucao manual (indices MATLAB)
caminhos = { [1 2 4 6], [1 3 5 6], [1 3 5 4 6] };   % 0-1-3-5, 0-2-4-5, 0-2-4-3-5
R = C;  total = 0;
fprintf('--- Caminhos da resolucao manual ---\n');
for k = 1:numel(caminhos)
    p = caminhos{k};
    g = inf;
    for i = 1:numel(p)-1, g = min(g, R(p(i), p(i+1))); end
    for i = 1:numel(p)-1
        R(p(i), p(i+1)) = R(p(i), p(i+1)) - g;
        R(p(i+1), p(i)) = R(p(i+1), p(i)) + g;
    end
    total = total + g;
    fprintf('P%d: %-16s gargalo = %2g | fluxo acumulado = %g\n', k, strjoin(nomes(p),'->'), g, total);
end
alcance = bfs_residual(R, s);
fprintf('Ainda existe caminho aumentante? %s\n', simnao(alcance(t)));
fprintf('Corte minimo: S = {%s} | T = {%s}\n', strjoin(nomes(alcance),','), strjoin(nomes(~alcance),','));

Fman = max(C - R, 0);
fprintf('\n%-8s %-6s %-6s %s\n','Aresta','c','f','c_f');
[I, J] = find(C');                       % percorre por linha (origem)
for e = 1:numel(I)
    u = J(e);  v = I(e);
    fprintf('%-8s %-6g %-6g %g\n', [nomes{u} '->' nomes{v}], C(u,v), Fman(u,v), C(u,v)-Fman(u,v));
end

%% 2) Implementacao do material (DFS; caminhos podem diferir, valor final nao)
[fluxoMax, F] = ford_fulkerson(C, s, t);
fprintf('\nFluxo maximo (ford_fulkerson do material): %g\n', fluxoMax);

%% 3) Funcao nativa
G = digraph(C, nomes);
[mf, GF, cs, ct] = maxflow(G, '0', '5');
fprintf('Fluxo maximo (maxflow nativo): %g\n', mf);
fprintf('Corte (maxflow nativo): S = {%s} | T = {%s}\n',strjoin(cs(:)', ','), strjoin(ct(:)', ','));
fprintf('Confere com a resolucao manual (23): %s\n', simnao(total == 23 && fluxoMax == 23 && mf == 23));

%% Figura: capacidades + fluxo da resolucao manual
figure('Name','Q4 - Ford-Fulkerson','Color','w');
rot = arrayfun(@(e) sprintf('%g/%g', Fman(findnode(G,G.Edges.EndNodes{e,1}), findnode(G,G.Edges.EndNodes{e,2})), ...
      G.Edges.Weight(e)), 1:numedges(G), 'UniformOutput', false);
Pl = plot(G,'EdgeLabel',rot,'NodeFontSize',14,'NodeColor','k','EdgeColor','b', ...
    'LineWidth',1.5,'Marker','o','MarkerSize',8,'ArrowSize',12,'Layout','layered','Sources',1,'Sinks',6);
satur = find(arrayfun(@(e) Fman(findnode(G,G.Edges.EndNodes{e,1}), findnode(G,G.Edges.EndNodes{e,2})) ...
      == G.Edges.Weight(e), 1:numedges(G)));
highlight(Pl, 'Edges', satur, 'EdgeColor', 'r', 'LineWidth', 4);
highlight(Pl, cs, 'NodeColor', 'g');
title('Q4 - Fluxo/capacidade (vermelho = saturada, verde = lado S do corte)');
salvar_fig(gcf, '../figs/q4_ford_fulkerson.png');
diary off;

%% ===================== FUNCOES =====================
function [fluxoMax, F] = ford_fulkerson(C, s, t)
    % Fluxo maximo: caminhos aumentantes por busca em profundidade no residual
    n = size(C, 1);  R = C;  fluxoMax = 0;
    while true
        pai = zeros(1, n);  pai(s) = s;  pilha = s;
        while ~isempty(pilha) && pai(t) == 0
            u = pilha(end);  pilha(end) = [];
            for v = find(R(u, :) > 0)
                if pai(v) == 0
                    pai(v) = u;  pilha(end+1) = v;
                end
            end
        end
        if pai(t) == 0, break; end
        gargalo = inf;  v = t;
        while v ~= s
            u = pai(v);  gargalo = min(gargalo, R(u, v));  v = u;
        end
        v = t;
        while v ~= s
            u = pai(v);
            R(u, v) = R(u, v) - gargalo;
            R(v, u) = R(v, u) + gargalo;
            v = u;
        end
        fluxoMax = fluxoMax + gargalo;
    end
    F = max(C - R, 0);
end

function vis = bfs_residual(R, s)
    % bfs() do material aplicada ao grafo residual (arestas com R > 0)
    n = size(R, 1);  vis = false(1, n);
    fila = zeros(1, n);  cab = 1;  cau = 1;
    fila(1) = s;  vis(s) = true;
    while cab <= cau
        v = fila(cab);  cab = cab + 1;
        for w = find(R(v, :) > 0)
            if ~vis(w)
                vis(w) = true;  cau = cau + 1;  fila(cau) = w;
            end
        end
    end
end

function s = simnao(x)
    if x, s = 'SIM'; else, s = 'NAO'; end
end

function salvar_fig(f, arq)
    if exist('exportgraphics','file'), exportgraphics(f, arq, 'Resolution', 200);
    else, saveas(f, arq); end
end