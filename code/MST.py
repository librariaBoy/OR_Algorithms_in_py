# 代码 6：图与网络优化 —— 最小生成树
import numpy as np

def main():
    '''
    定义问题形式：
    1. 节点数: n
    2. 边数: m
    3. 边集合: edges
    4. 边的权值: weights
    5. 邻接矩阵: A
    '''
    n = 5
    m = 8

    A = np.array([ [0, 1, 2, 0, 0],
                   [1, 0, 2, 4, 3],
                   [2, 2, 0, 4, 4],
                   [0, 4, 4, 0, 2],
                   [0, 3, 4, 2, 0]])

    edges = {(0,1):1,
             (0,2):2,
             (1,2):2,
             (1,3):4,
             (1,4):3,
             (2,3):4,
             (2,4):4,
             (3,4):2}
    
    MST = Graph(n, m, edges, A)
    S_1 = MST.Kruskal()
    print("Kruskal:")
    print(S_1)

    print("")
    S_2 = MST.Dijkstra()
    print("Dijkstra:")
    print(S_2)    
    


class Graph:
    '''
    图类
    '''

    def __init__(self, n, m, edges, A):
        '''
        初始化图
        :param n: 节点数
        :param m: 边数
        :param A: 邻接矩阵
        '''
        self.n = n
        self.m = m
        self.edges = edges
        self.A = A

        # 检查是否成环
        self.Occupied = [i for i in range(self.n)]

    def Kruskal(self):
        '''
        Kruskal算法
        :return: 最小生成树
        '''
        # 从小到大整理边
        sort_edges = dict(sorted(self.edges.items(), key=lambda x: x[1]))

        # 存放最小树的边
        S = []
        i = 0
        j = 1

        for edge in sort_edges:
            # 结束条件
            if len(S) == self.n - 1:
                return S
            
            # 检查是否构成回路
            if self.not_ring(edge):
                # 如果没成环，改变标号，并添加此边
                self.change_mark(edge,S)
                S.append(edge)
                

                i = i + 1
                j = j + 1
        
        j = j + 1

    def Dijkstra(self):
        '''
        Dijkstra算法
        :return: 最小生成树
        '''

        T = []
        R = {0}
        S = {i for i in range(self.n) if i not in R}

        mark = [[i, float('inf')] for i in range(self.n)]

        while len(S) != 0:

            mark = self.gen_mark(R, S, mark)

            # 找的 mark 中最小的 加入 R
            index = mark.index(min(mark, key=lambda x: x[1]))

            R.add(index)
            mark[index][1] = float('inf')
            T.append((mark[index][0], index))
            S.remove(index)

        return T
                    

            

    def not_ring(self, edge):
        '''
        检查是否形成回路
        '''
        if self.Occupied[edge[0]] != self.Occupied[edge[1]]:
            return True
        else:
            return False
        
    def change_mark(self, edge, S):
        '''
        改变标号
        '''
        flag = False
        for cur_edge in S:
            if cur_edge[0] in edge or cur_edge[1] in edge:
                flag = True
                break

        if flag:
            # 若已加入，之前的符号全部改变
            self.Occupied[edge[0]] = self.Occupied[edge[1]]
            for cur_edge in S:
                self.Occupied[cur_edge[0]] = self.Occupied[edge[0]]
                self.Occupied[cur_edge[1]] = self.Occupied[edge[0]]
        else:
            # 若未加入,只改变当前段的
            self.Occupied[edge[0]] = self.Occupied[edge[0]]
            
    def gen_mark(self, R, S, mark):
        '''
        生成标号
        '''
        for i in R:
            for j in S:

                # 如果边存在
                if (i,j) in self.edges:
                    if mark[j][1] > self.edges[(i,j)]:
                        mark[j][0] = i
                        mark[j][1] = self.edges[(i,j)]
                elif (j,i) in self.edges:
                    if mark[j][1] > self.edges[(j,i)]:
                        mark[j][0] = i
                        mark[j][1] = self.edges[(j,i)]
        return mark
            
main()       