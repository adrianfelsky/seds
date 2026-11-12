% q3_prim.m - Validacao da 3a Questao (Prim, inicio em A)
% Algoritmo: prim() do material "algoritmos em MATLAB para Grafos"
%            + minspantree nativo (Method 'dense' = Prim), como no Lab. 01.
% Requer MATLAB R2016b ou superior.
clear; close all; clc;
if ~exist('../figs','dir'), mkdir('../figs'); end
if ~exist('../out','dir'),  mkdir('../out');  end
arqOut = '../out/q3_prim.txt';
if exist(arqOut,'file'), delete(arqOut); end
diary(arqOut);

%% Grafo nao direcionado (matriz simetrica)
nomes = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
L = 'ABCDEFG';
arestas = { 'A','D', 5;  'A','B', 7;  'B','C', 8;  'B','D', 9;
            'B','E', 7;  'C','E', 5;  'D','E',15;  'D','F', 6;
            'E','F', 8;  'E','G', 9;  'F','G',11 };
P = zeros(7);
for e = 1:size(arestas,1)
    P(L == arestas{e,1}, L == arestas{e,2}) = arestas{e,3};
end
P = P + P';

%% Prim (implementacao do material, com impressao dos passos)
[ag, pesoTotal] = prim_tab(P, 1, nomes);
fprintf('\n--- AGM (ordem de inclusao) ---\n');
for k = 1:size(ag,1)
    fprintf('%da: (%s, %s) = %g\n', k, nomes{ag(k,1)}, nomes{ag(k,2)}, P(ag(k,1),ag(k,2)));
end
fprintf('Peso total: %g\n', pesoTotal);

%% Funcao nativa
G = graph(P, nomes);
[T, ~] = minspantree(G, 'Method', 'dense', 'Root', 'A');
fprintf('\nPeso total (minspantree nativo): %g\n', sum(T.Edges.Weight));
fprintf('Confere com a resolucao manual (39): %s\n', simnao(pesoTotal == 39 && sum(T.Edges.Weight) == 39));

%% Figura
figure('Name','Q3 - Prim','Color','w');
Pl = plot(G,'EdgeLabel',G.Edges.Weight,'NodeFontSize',14,'NodeColor','k', ...
    'EdgeColor','b','LineWidth',1.5,'Marker','o','MarkerSize',8,'Layout','force');
highlight(Pl, T, 'EdgeColor','g', 'LineWidth', 4);
title('Q3 - Arvore Geradora Minima (Prim, inicio A)');
salvar_fig(gcf, '../figs/q3_prim.png');
diary off;

%% ===================== FUNCOES =====================
function [arestas, pesoTotal] = prim_tab(W, inicio, nomes)
    % prim() do material + impressao de chave/pai a cada passo
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
            else, txt = sprintf('%g(%s)', custo(v), nomes{pai(v)}); end
            fprintf('%-8s', txt);
        end
        fprintf('\n');
    end
end

function s = simnao(x)
    if x, s = 'SIM'; else, s = 'NAO'; end
end

function salvar_fig(f, arq)
    if exist('exportgraphics','file'), exportgraphics(f, arq, 'Resolution', 200);
    else, saveas(f, arq); end
end