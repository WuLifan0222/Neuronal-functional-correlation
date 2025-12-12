function exists = hasEdge(G, node1, node2)
% 检查两个节点之间是否有边
exists = any(ismember(G.Edges.EndNodes, [node1, node2], 'rows'));
