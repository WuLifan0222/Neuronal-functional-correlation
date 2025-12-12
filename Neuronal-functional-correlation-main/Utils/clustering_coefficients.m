function C = clustering_coefficients(G)
n = numnodes(G); % 节点数
C = zeros(1, n); % 初始化聚集系数数组
for i = 1:n
    nbs = neighbors(G, i); % 获取节点 i 的邻居
    k = length(nbs); % 邻居的数量
    if k < 2
        C(i) = 0; % 如果邻居少于2，聚集系数为0
    else
        edgesInTriangle = 0; % 三角形中的边数
        for j = 1:k
            for z = j+1:k
                if any(ismember(G.Edges.EndNodes, [nbs(j), nbs(z)], 'rows'))
                    edgesInTriangle = edgesInTriangle + 1;
                end
            end
        end
        C(i) = (2 * edgesInTriangle) / (k * (k - 1)); % 计算聚集系数
    end
end
