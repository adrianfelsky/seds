% run_all.m - Executa todas as validacoes da Prova 2 (Grafos).
% Rode a partir da pasta "matlab/". Gera ../figs/*.png e ../out/*.txt,
% que sao incluidos automaticamente no documento LaTeX.
q1_dijkstra;
q2_bellman_ford;
q3_prim;
q4_ford_fulkerson;
q5_kruskal;
disp('Validacoes concluidas: figuras em ../figs e saidas em ../out');