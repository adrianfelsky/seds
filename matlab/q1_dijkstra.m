clear; close all; clc;
if ~exist('../figs','dir'), mkdir('../figs'); end
if ~exist('../out','dir'),  mkdir('../out');  end
arqOut = '../out/q1_dijkstra.txt';
if exist(arqOut,'file'), delete(arqOut); end
diary(arqOut);

nomes = {'0' '1' '2' '3' '4' '5'};
W = zeros(6);
W(1,3) = 20;  W(1,4) = 40;  W(1,5) = 30;   % 0->2, 0->3, 0->4
W(2,3) = 10;                               % 1->2
W(3,5) = 15;                               % 2->4
W(4,5) = 10;  W(4,6) = 80;                 % 3->4, 3->5
W(5,2) = 82;  W(5,6) = 20;                 % 4->1, 4->5

%% Dijkstra
[dist, pred] = dijkstra_tab(W, 1, nomes);

fprintf('\n--- Resultado final (origem 0) ---\n');
fprintf('%-8s %-8s %-12s %s\n','Destino','Custo','Antecessor','Caminho');
for v = 1:6
    p = caminho(pred, 1, v);
    if pred(v) == 0, a = '--'; else, a = nomes{pred(v)}; end
    fprintf('%-8s %-8g %-12s %s\n', nomes{v}, dist(v), a, strjoin(nomes(p),' -> '));
end

%% Validacao com funcoes nativas do MATLAB
G = digraph(W, nomes);
dNat = distances(G, '0');
fprintf('\nDistancias (distances nativo): %s\n', mat2str(dNat));
esperado = [0 112 20 40 30 50];                    % resolucao manual
fprintf('Confere com a resolucao manual: %s\n', simnao(isequal(dist, esperado) && isequal(dNat, esperado)));

%% Figura
figure('Name','Q1 - Dijkstra','Color','w');
Pl = plot(G,'EdgeLabel',G.Edges.Weight,'NodeFontSize',14,'NodeColor','k', ...
    'EdgeColor','b','LineWidth',1.5,'Marker','o','MarkerSize',8,'ArrowSize',12,'Layout','layered');
TR = shortestpathtree(G, '0');
highlight(Pl, TR, 'EdgeColor','g', 'LineWidth', 4);
title('Q1 - Arvore de caminhos minimos (Dijkstra, origem 0)');
salvar_fig(gcf, '../figs/q1_dijkstra.png');
diary off;

%% ===================== FUNCOES =====================
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
            else, txt = sprintf('%g(%s)', dist(v), nomes{pred(v)}); end
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

function s = simnao(x)
    if x, s = 'SIM'; else, s = 'NAO'; end
end

function salvar_fig(f, arq)
    if exist('exportgraphics','file'), exportgraphics(f, arq, 'Resolution', 200);
    else, saveas(f, arq); end
end