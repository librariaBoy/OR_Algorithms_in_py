# 代码 1 线性规划
import numpy as np

def main():
    '''
    在此编辑问题形式
    min z = c^T x
    s.t. Ax = b
         x >= 0
    '''
    c = np.array([1., 1., 1., 0., 0.])

    A = np.array([[-3., -1., -1., 1., 0.],
                  [ 1., -4., -1., 0., 1.]])
    
    b = np.array([-1., -2.])
    
    '''调用单纯形法'''
    # [x,obj_value] = simplex_solver(A, b, c)

    '''调用对偶单纯形法'''
    [x,obj_value] = dual_simplex_solver(A, b, c)

    print("目标函数值:",obj_value)
    print("最优解:",x)

def simplex_solver(A, b, c):
    '''
    单纯形算法
    必须在A中有基本可行基
    A: 约束矩阵
    b: 右端向量
    c: 价值向量
    '''

    print("simplex_solver v_1.0")

    # 预处理
    # 识别该问题的一些参数

    # m 行 n 列
    m, n = A.shape

    # 基变量
    base = list(range(n-m, n))

    # 检验数向量
    zeta = -c

    # 目标值
    obj_value = 0

    # 转变为初始单纯性表
    [zeta,obj_value] = trans_sheet(A,b,zeta,obj_value,base)

    # 主循环
    itr = 0     #迭代次数
    
    while True:
        # 迭代次数提示
        itr = itr + 1
        print("iteration:",itr)

        # step 3 求最大的检验数
        max_zeta = max(zeta)

        ## step 4 判断是否满足最优
        if max_zeta <= 0:
            # 根据基变量得出最优解
            x = []
            for i in range(n):
                if i in base:
                    x.append(b[base.index(i)])
                else:
                    x.append(0)

            return [x,obj_value]
            break

        ## step 5 求出最大的检验数对应的列标号, 检查是否该列小于0
        k = zeta.tolist().index(max_zeta)
        if (A[:,k] <= 0).all():
            print("无界")
            break

        ## step 6求出该列最小的 b/a， 确定出基位置(a > 0)
        r = 0   # 记录出基位置
        min_ratio = 10000000
        for i in range(m):
            if A[i,k] <= 0:
                continue

            # 只寻找 a > 0 的行
            if b[i]/A[i,k] < min_ratio:
                min_ratio = b[i]/A[i,k]
                r = i

        ## step 7 替换基变量
        base[r] = k
        # 化为标准单纯性表
        b[r] = b[r]/A[r,k]
        A[r,:] = A[r,:]/A[r,k]
        # 把这一列所在行消灭
        for i in range(m):
            if i != r:
                b[i] = b[i] - A[i,k]*b[r]
                A[i,:] = A[i,:] - A[i,k]*A[r,:]
        [zeta,obj_value] = trans_sheet(A,b,zeta,obj_value,base)


def dual_simplex_solver(A, b, c):
    '''
    对偶单纯形算法
    必须含有一个基本解和一个对偶问题的可行解
    A: 约束矩阵
    b: 右端向量
    c: 价值向量
    '''
    
    print("dual_simplex_solver v_1.0")

    # 预处理
    m, n = A.shape
    base = list(range(n-m, n))

    # 目标函数值
    obj_value = np.array([0])

    # 直接利用函数得到单纯性表
    zeta = np.hstack((-c,obj_value))
    A = np.hstack((A,b.reshape(-1,1))) 

    itr = 0
    while True:
        itr = itr + 1
        print("iteration:",itr)
        # step 2 求最小的 b 值
        b_r = A[:,-1].min()
        r = A[:,-1].tolist().index(b_r)

        # step 3 判断最优性
        if b_r >= 0:
            # 根据基变量得出最优解
            x = []
            for i in range(n):
                if i in base:
                    x.append(A[base.index(i),-1])
                else:
                    x.append(0)

            return [x,zeta[-1]]

        # step 4 判断无解
        if (A[r,:] >= 0).all():
            print("无解")
            break

        # step 5 求出该列最小的 zeta/a， 确定入基位置(a_rj < 0)
        k = 0   # 记录出基位置
        min_ratio = 10000000
        for j in range(n):
            if A[r,j] >= 0:
                continue

            # 只寻找 a > 0 的行
            if A[r,-1]/A[r,j] < min_ratio:
                min_ratio = A[r,-1]/A[r,j]
                k = j

        # step 6 替换基变量
        base[r] = k
        A[r,:] = A[r,:]/A[r,k]
        for i in range(m):
            if i != r:
                A[i,:] = A[i,:] - A[i,k]*A[r,:]
        zeta = zeta - A[r,:]*zeta[k]




def trans_sheet(A,b,zeta,obj_value,base):
    '''
    将存在基本可行基的约束矩阵A转换为初始单纯性表
    '''
    # 行标号
    i = 0
    # 行变换
    for j in base:
        obj_value = obj_value - zeta[j]*b[i]
        zeta = zeta - A[i,:]*zeta[j]
        i = i + 1
    
    return [zeta,obj_value]


main()