% q5_kruskal.m - Validacao da 5a Questao (Kruskal)
% Algoritmo: Kruskal com Union-Find conforme pseudocodigo do material (1.pdf)
%            + minspantree nativo (Method 'sparse' = Kruskal), como no Lab. 01.
% Rotulos: A (sup. esq.), B (sup. centro), C (sup. dir.),
%          D (inf. esq.), E (base), F (inf. dir.)
% Requer MATLAB R2016b ou superior.
clear; close all; clc;
if ~exist('../figs','dir'), mkdir('../figs'); end
if ~exist('../out','dir'),  mkdir('../out');  end
arqOut = '../out/q5_kruskal.txt';
if exist(arqOut,'file'), delete(arqOut); end
diary(arqOut);

%% Grafo nao direcionado
nomes = {'A' 'B' 'C' 'D' 'E' 'F'};
L = 'ABCDEF';
arestas = { 'A','B',2;  'B','C',3;  'A','D',5;  'A','E',1;  'B','E',7;
            'C','E',3;  'C','F',9;  'D','E',4;  'E','F',4 };
W = zeros(6);
for e = 1:size(arestas,1)
    W(L == arestas{e,1}, L == arestas{e,2}) = arestas{e,3};
end
W = W + W';

%% Kruskal (pseudocodigo do material, com impressao das decisoes)
[ag, pesoTotal] = kruskal_tab(W, nomes);
fprintf('\nArestas da AGM: ');
par = [nomes(ag(:,1)); nomes(ag(:,2))];
fprintf('(%s,%s) ', par{:});
fprintf('\nPeso total: %g\n', pesoTotal);

%% Funcao nativa
G = graph(W, nomes);
T = minspantree(G, 'Method', 'sparse');
fprintf('Peso total (minspantree sparse nativo): %g\n', sum(T.Edges.Weight));
fprintf('Confere com a resolucao manual (14): %s\n', simnao(pesoTotal == 14 && sum(T.Edges.Weight) == 14));

%% Figura (posicoes semelhantes ao desenho da prova)
X = [1.5 3.0 4.3 1.2 2.8 4.1];
Y = [3.0 3.2 3.0 1.5 1.0 1.5];
figure('Name','Q5 - Kruskal','Color','w');
Pl = plot(G,'XData',X,'YData',Y,'EdgeLabel',G.Edges.Weight,'NodeFontSize',14, ...
    'NodeColor','k','EdgeColor','b','LineWidth',1.5,'Marker','o','MarkerSize',8);
highlight(Pl, T, 'EdgeColor','g', 'LineWidth', 4);
title('Q5 - Arvore Geradora Minima (Kruskal)');
salvar_fig(gcf, '../figs/q5_kruskal.png');
diary off;

%% ===================== FUNCOES =====================
function [A, pesoTotal] = kruskal_tab(W, nomes)
    % Kruskal(G): ordena arestas por peso e usa Union-Find para evitar ciclos
    n = size(W, 1);
    [u, v] = find(triu(W));
    w = W(sub2ind(size(W), u, v));
    E = sortrows([w u v], [1 2 3]);       % peso crescente (empate: ordem alfabetica)
    pai = 1:n;                            % Criar-Conjunto(v)
    A = zeros(0, 2);  pesoTotal = 0;
    fprintf('Passo | Aresta | Peso | Decisao\n');
    for k = 1:size(E, 1)
        a = E(k,2);  b = E(k,3);
        ra = buscar(pai, a);  rb = buscar(pai, b);
        if ra ~= rb
            A(end+1, :) = [a b];  pesoTotal = pesoTotal + E(k,1);
            pai(rb) = ra;                 % Unir-Conjuntos(u, v)
            dec = 'Aceita';
        else
            dec = 'Rejeitada (ciclo)';
        end
        fprintf('%5d | (%s,%s)  | %4g | %s\n', k, nomes{a}, nomes{b}, E(k,1), dec);
    end
end

function r = buscar(pai, x)
    % Buscar-Conjunto(x)
    r = x;
    while pai(r) ~= r, r = pai(r); end
end

function s = simnao(x)
    if x, s = 'SIM'; else, s = 'NAO'; end
end

function salvar_fig(f, arq)
    if exist('exportgraphics','file'), exportgraphics(f, arq, 'Resolution', 200);
    else, saveas(f, arq); end
end