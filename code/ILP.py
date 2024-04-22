# 代码 2 整数规划
import math
import numpy as np
import gurobipy as gp

def main():
    '''
    在此编辑问题形式
    min z = c^T x
    s.t. Ax = b
         x >= 0 且为整数
    '''

    # 松弛前问题决策变量数
    num = 2

    c = np.array([-1., -1., 0., 0., 0.])

    A = np.array([[-4., 2., 1., 0., 0.],
                  [ 4., 2., 0., 1., 0.],
                  [ 0.,-2., 0., 0., 1.]])
    
    b = np.array([-1., 11., -1.])
    
    '''调用分支定界法'''
    root = branch_and_bound_solver(A, b, c)

    print("                         ")
    print("-------------------------")
    print("                         ")
    print("目标函数值:",root.obj_value)
    print("最优解:",root.x[:num])
    print("                         ")
    print("-------------------------")
    print("                         ")
    
# 二叉树的节点类
class Node:
    def __init__(self, A, b, c):

        # 节点对应问题
        self.A = A
        self.b = b
        self.c = c

        # 节点对应解
        self.obj_value = float('inf')
        self.x = None

    def solve(self):
        '''
        调用 Gurobi 求解松弛问题
        '''
        SP = gp.Model('SP')
        SP.setParam('OutputFlag', 0)
        x = SP.addVars(len(self.c), lb=0, vtype=gp.GRB.CONTINUOUS)
        SP.setObjective(gp.quicksum(self.c[i]*x[i] for i in range(len(self.c))), gp.GRB.MINIMIZE)
        SP.addConstrs((gp.quicksum(self.A[i][j]*x[j] for j in range(len(self.A[i]))) == self.b[i] for i in range(len(self.b))))
        SP.optimize()

        if SP.status == gp.GRB.Status.OPTIMAL:
            self.x = [x[i].x for i in range(len(x))]
            self.obj_value = SP.objVal

            if self.obj_value < self.current_Obj:
                print("        当前最优值为:", self.obj_value)
                print("        当前最优解为:", self.x)
            else:
                self.cut = False

        else:
            print("        无解")
                

    def check_not_integer(self):
        '''
        返回解不是整数的第一项的索引
        '''
        for i in range(len(self.x)):
            if self.x[i] % 1 != 0:
                return i
        return -1


def branch_and_bound_solver(A, b, c):
    '''
    分支定界法
    A: 矩阵
    b: 向量
    c: 向量
    '''

    print("branch_and_bound_solver v_1.0")
    print("------------------")

    print("根节点")

    root = Node(A, b, c)
    root.solve()

    # 当前最优值
    current_Obj = root.obj_value

    
    i = 1
    while True:
        # 分支
        index = root.check_not_integer()

        if index == -1:
            return root

        print("------------------")
        print("分支树 第", i, "层")

        # 左子树操作
        print("    左子树:")
        [new_A, new_b, new_c] = gen_matrix(root.A, root.b, root.c, index, root.x[index], True)
        root_left = Node(new_A, new_b, new_c)
        root_left.solve()

        # 右子树操作
        print("    右子树:")
        [new_A, new_b, new_c] = gen_matrix(root.A, root.b, root.c, index, root.x[index], False)
        root_right = Node(new_A, new_b, new_c)
        root_right.solve()

        if root_left.obj_value < root_right.obj_value:
            root = root_left
        else:
            root = root_right

        i = i + 1



def gen_matrix(A, b, c, index, key, flag):
    '''
    生成添加新约束后的矩阵
    key: 最优解的值
    flag = true 为  <=
    flag = false 为 >=
    '''
    m,n = A.shape
    new_A_col = np.zeros((m+1,1))
    new_A_col[m,0] = 1


    new_A_row = np.zeros((1,n))
    if flag:
        new_A_row[0,index] = 1
        b = np.append(b,math.floor(key))
    else:
        new_A_row[0,index] = -1
        b = np.append(b, -math.ceil(key))

    A = np.vstack((A, new_A_row))
    A = np.hstack((A, new_A_col))

    c = np.append(c, 0)

    return A, b, c


main()