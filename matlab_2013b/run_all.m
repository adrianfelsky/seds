% run_all.m - Executa todas as validacoes da Prova 2 (Grafos). MATLAB R2013b.
% Rode com a pasta "matlab/" como Current Folder.
% Gera ../figs/*.png e ../out/*.txt, incluidos automaticamente no LaTeX.
q1_dijkstra;
q2_bellman_ford;
q3_prim;
q4_ford_fulkerson;
q5_kruskal;
disp('Validacoes concluidas: figuras em ../figs e saidas em ../out');