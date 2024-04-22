# 代码 7：图与网络优化 —— 最短有向路

def main():
    '''
    定义问题形式：
    1. 节点数: n
    2. 边数: m
    3. 边集合: edges
    4. 边的权值: weights
    5. 邻接矩阵: A
    '''
    n = 6

    edges = {(0,1):5,
             (0,3):3,
             (0,5):4,
             (1,2):3,
             (1,3):2,
             (2,5):2,
             (2,4):4,
             (3,2):7,
             (3,4):5,
             (4,1):4,
             (4,2):1,
             (5,4):6}
    
    # 设置起点
    start = 1

    
    graph = Graph(n, edges, start)

    print("")
    distance, path = graph.Dijkstra()
    # 结果展示
    print("Dijkstra:")
    print("distance=",distance)

    for i in range(n):
        if distance[i] == float("inf"):
            print("no path from",start,"to",i)
        else:
            print("path from",start,"to",i,"is:")
            for j in path[i]:
                print(j,"-> ", end="")
            print(i)



class Graph:
    '''
    图类
    '''

    def __init__(self, n, edges, start):
        '''
        初始化图
        :param n: 节点数
        :param edges: 边
        :param start: 起点
        '''
        self.n = n
        self.edges = edges
        self.start = start


    def Dijkstra(self):
        '''
        Dijkstra算法
        :return: 最短有向路
        '''
        # 路径记录
        path = [[self.start] for _ in range(self.n)]
        
        ## step 1 开始
        # 初始化
        u = [0]*self.n
        P = {self.start}
        T = {i for i in range(self.n) if i not in P}

        # 生成初始 u 列表
        itr = 1
        for i in range(self.n):
            if i == self.start:
                u[i] = 0
            else:
                if (self.start, i) in self.edges:
                    u[i] = self.edges[(self.start, i)]
                else:
                    u[i] = float('inf')

        while True:

            if T == set():
                return u, path
            
            print("第", itr, "轮循环")
            
            # step 2 永久标号         
            k = min(T, key=lambda i: u[i])

            P.add(k)
            T.remove(k)
            

            # step 3 修改临时标号
            for j in range(self.n):
                if (k,j) in self.edges:
                    if u[j] > u[k] + self.edges[(k,j)]:
                        u[j] = u[k] + self.edges[(k,j)]
                        path[j].append(k)
            
            itr = itr + 1
            
main()       