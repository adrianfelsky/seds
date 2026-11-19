function f = desenha_grafo(X, Y, nomes, W, dirigido, D, rot, titulo, corNo)
% DESENHA_GRAFO  Desenha um grafo sem usar graph/digraph (compativel R2013b).
%   X, Y     : coordenadas dos vertices
%   nomes    : cell com os rotulos dos vertices
%   W        : matriz de pesos/capacidades (0 = sem aresta)
%   dirigido : true (setas) ou false
%   D        : matriz logica de arestas destacadas (verde e grossa)
%   rot      : cell n x n com rotulo de cada aresta ('' = usa o peso)
%   titulo   : titulo da figura
%   corNo    : matriz n x 3 com a cor de cada vertice
n = numel(X);
if nargin < 6 || isempty(D),     D = false(n);  end
if nargin < 7 || isempty(rot),   rot = cell(n); end
if nargin < 8,                   titulo = '';   end
if nargin < 9 || isempty(corNo), corNo = repmat([0.75 0.85 1], n, 1); end

esc = max(max(X) - min(X), max(Y) - min(Y));
if esc == 0, esc = 1; end
r = 0.06 * esc;                         % raio do vertice

f = figure('Color', 'w');
hold on; axis equal; axis off;

if dirigido
    M = W;  DD = D;
else
    M = triu(W);  DD = triu(D | D');
end
[I, J] = find(M);

for e = 1:numel(I)
    i = I(e);  j = J(e);
    p1 = [X(i) Y(i)];  p2 = [X(j) Y(j)];
    d = p2 - p1;  u = d / norm(d);  nrm = [-u(2) u(1)];
    off = 0;
    if dirigido && W(j, i) ~= 0, off = 0.05 * esc; end   % par i->j e j->i
    a = p1 + u * r + nrm * off;
    b = p2 - u * r + nrm * off;
    if DD(i, j), cor = [0 0.7 0]; lw = 4; else cor = [0 0 1]; lw = 1.5; end
    plot([a(1) b(1)], [a(2) b(2)], '-', 'Color', cor, 'LineWidth', lw);
    if dirigido
        ha = 0.07 * esc;  wa = 0.03 * esc;
        base = b - u * ha;
        pts = [b; base + nrm * wa; base - nrm * wa];
        patch('XData', pts(:,1), 'YData', pts(:,2), 'FaceColor', cor, 'EdgeColor', cor);
    end
    m = (a + b) / 2 + nrm * 0.045 * esc;
    if isempty(rot{i, j}), txt = num2str(W(i, j)); else txt = rot{i, j}; end
    text(m(1), m(2), txt, 'HorizontalAlignment', 'center', 'FontSize', 10, ...
        'BackgroundColor', 'w', 'Margin', 1, 'Interpreter', 'none');
end

for k = 1:n
    rectangle('Position', [X(k)-r, Y(k)-r, 2*r, 2*r], 'Curvature', [1 1], ...
        'FaceColor', corNo(k, :), 'EdgeColor', 'k', 'LineWidth', 1.5);
    text(X(k), Y(k), nomes{k}, 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');
end
xlim([min(X) - 2*r, max(X) + 2*r]);
ylim([min(Y) - 2*r, max(Y) + 2*r]);
title(titulo, 'Interpreter', 'none');
hold off;
end