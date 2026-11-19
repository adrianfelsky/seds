function q4_ford_fulkerson
% Q4_FORD_FULKERSON - Validacao da 4a Questao (s = 0, t = 5). MATLAB R2013b.
% Algoritmo: ford_fulkerson() do material + caminhos da resolucao manual.
% Conferencia extra: graphmaxflow (Bioinformatics Toolbox).
close all; clc;
prepara_pastas();
arqOut = '../out/q4_ford_fulkerson.txt';
if exist(arqOut, 'file'), delete(arqOut); end
diary(arqOut);

% Rede (vertice k da prova = indice k+1 no MATLAB)
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

% 1) Caminhos aumentantes da resolucao manual
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
    fprintf('P%d: %-16s gargalo = %2g | fluxo acumulado = %g\n', ...
        k, strjoin(nomes(p), '->'), g, total);
end
alcance = bfs_residual(R, s);
fprintf('Ainda existe caminho aumentante? %s\n', simnao(alcance(t)));
fprintf('Corte minimo: S = {%s} | T = {%s}\n', ...
    strjoin(nomes(alcance), ','), strjoin(nomes(~alcance), ','));

Fman = max(C - R, 0);
fprintf('\n%-8s %-6s %-6s %s\n', 'Aresta', 'c', 'f', 'c_f');
for u = 1:6
    for v = find(C(u, :))
        fprintf('%-8s %-6g %-6g %g\n', [nomes{u} '->' nomes{v}], C(u,v), Fman(u,v), C(u,v) - Fman(u,v));
    end
end

% 2) Implementacao do material (DFS; caminhos podem diferir, valor final nao)
fluxoMax = ford_fulkerson(C, s, t);
fprintf('\nFluxo maximo (ford_fulkerson do material): %g\n', fluxoMax);

ok = (total == 23) && (fluxoMax == 23);
if tem_bioinfo()
    [mf, FM, corte] = graphmaxflow(sparse(C), s, t); %#ok<ASGLU>
    ladoS = corte(1, :) == 1;
    fprintf('Fluxo maximo (graphmaxflow): %g\n', mf);
    fprintf('Corte (graphmaxflow): S = {%s} | T = {%s}\n', ...
        strjoin(nomes(ladoS), ','), strjoin(nomes(~ladoS), ','));
    ok = ok && (mf == 23);
else
    fprintf('(Bioinformatics Toolbox indisponivel: conferencia so com a resolucao manual)\n');
end
fprintf('Confere com a resolucao manual (23): %s\n', simnao(ok));

% Figura: rotulo fluxo/capacidade, saturadas em destaque, lado S em verde
rot = cell(6);
for u = 1:6
    for v = find(C(u, :))
        rot{u, v} = sprintf('%g/%g', Fman(u, v), C(u, v));
    end
end
D = (C > 0) & (Fman == C);
corNo = repmat([0.75 0.85 1], 6, 1);
corNo(alcance, :) = repmat([0.6 1 0.6], sum(alcance), 1);
X = [0 1 1 2  2 3];
Y = [0 1 -1 1 -1 0];
f = desenha_grafo(X, Y, nomes, C, true, D, rot, ...
    'Q4 - Fluxo/capacidade (verde grosso = saturada; vertice verde = lado S)', corNo);
print(f, '-dpng', '-r200', '../figs/q4_ford_fulkerson.png');
diary off;
end

% ===================== SUBFUNCOES =====================
function fluxoMax = ford_fulkerson(C, s, t)
n = size(C, 1);  R = C;  fluxoMax = 0;
while true
    pai = zeros(1, n);  pai(s) = s;  pilha = s;
    while ~isempty(pilha) && pai(t) == 0
        u = pilha(end);  pilha(end) = [];
        for v = find(R(u, :) > 0)
            if pai(v) == 0
                pai(v) = u;  pilha(end+1) = v; %#ok<AGROW>
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
end

function vis = bfs_residual(R, s)
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

function ok = tem_bioinfo()
ok = license('test', 'Bioinformatics_Toolbox') && exist('graphmaxflow', 'file') > 0;
end

function s = simnao(x)
if x, s = 'SIM'; else s = 'NAO'; end
end

function prepara_pastas()
if ~exist('../figs', 'dir'), mkdir('../figs'); end
if ~exist('../out', 'dir'),  mkdir('../out');  end
end