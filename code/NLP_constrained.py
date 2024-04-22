# 代码 5 非线性规划——含约束
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import LineSearch

def main():
    '''
    完成一些问题的设置
    和最后解的显示
    '''

    # 设置初始值
    x0 = np.array([0., 0., 2., 5.])

    # 设置约束形式
    A = np.array([[ 1., 1., 1., 0.], 
                  [ 1., 5., 0., 1.]])
    
    b = np.array([2., 5.])

    # 设置精度
    epsilon = 1e-6

    print("Wolfe_Solver v_1.0")

    # 调用 wolfe 法
    x_opt, x, y = wolfe_method(Objective, A, b, x0, epsilon)

    value = np.zeros(x.shape)

    for i in range(x.shape[0]):
        value[i] = Objective(np.array([x[i], y[i]]))

    print("x_opt=", x_opt[:b.shape[0]])
    print('f(x_opt) =', Objective(x_opt))

    # 可视化
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot(x, y, value, label='f(x)', color='blue', marker='x')
    
    plt.legend(loc = 'best')
    ax.set_xlabel('x_1')
    ax.set_ylabel('x_2')
    ax.set_zlabel('f(x)')

    plt.show()


def Objective(x):
    '''
    定义目标函数
    '''
    return 2*x[0]**2 + 2*x[1]**2 - 2*x[0]*x[1] - 4*x[0] - 6*x[1]

def wolfe_method(fun, A, b, x0, epsilon=1e-5):
    '''
    Wolfe 简约梯度法
    输入：
        fun: 目标函数
        A: 约束矩阵
        b: 右端向量
        x0: 初始值
        epsilon: 精度
    输出：
        x_opt: 最优解
        x: 迭代过程的 x 值
        y: 迭代过程的 f(x) 值
    '''
    # 预处理
    m,n = A.shape

    # step 1
    x = x0
    x_list = np.array(x[0])
    y_list = np.array(x[1])
    k = 0

    while True:
        # step 2
        # 最大分量下标集
        I_B, I_ext = gen_I_B(x,m)
        B_k, N_k = decomposition(A, I_B)

        # step 3
        grad = gradient(fun, x)
        grad_B, grad_N = grad_decomposition(grad, I_B)
        # 计算简约梯度
        r_N = -(np.linalg.inv(B_k) @ N_k).T @ grad_B + grad_N

        # step 4
        # 可行方向
        p = construct_p(B_k, N_k, r_N, I_B, I_ext, x)

        if norm2(p) <= epsilon:
            break
        else:
            dynamic_fun = set_dynamic_fun(x, p)
            ub = get_upper_bound(p, x)
            a,b = LineSearch.get_search_bound(dynamic_fun)
            t = LineSearch.golden_section_search(0,min([b,ub]),dynamic_fun, epsilon)
            x = x + t*p

            x_list = np.append(x_list, x[0])
            y_list = np.append(y_list, x[1])

            k = k + 1
            print("第",k, "次迭代")
    
    return x, x_list, y_list

def gen_I_B(x,m):
    '''
    根据 x, 生成 m 个最大分量下标集
    '''
    n = x.shape[0]
    I_B = np.argsort(-x)[:m].tolist()
    I_B.sort()
    I_ext = [i for i in range(n) if i not in I_B]

    return I_B, I_ext

def decomposition(A, I_B):
    '''
    根据I_B, 分解 A
    输入：
        A: 约束矩阵
        b: 右端向量
        I_B: 最大分量下标集
    输出：
        B_k: 矩阵 B_k
        N_k: 矩阵 N_k
    '''
    n = A.shape[1]
    B_k = A[:, I_B]
    N_k = A[:, [i for i in range(n) if i not in I_B]]
    return B_k, N_k

def grad_decomposition(grad, I_B):
    '''
    根据 I_B, 分解 grad
    输入：
        grad: 梯度
        I_B: 最大分量下标集
    输出：
        grad_B: 梯度 grad_B
        grad_N: 梯度 grad_N
    '''
    n = grad.shape[0]
    grad_B = grad[I_B]
    grad_N = grad[[i for i in range(n) if i not in I_B]]
    return grad_B, grad_N

def get_upper_bound(p, x):
    '''
    根据 p, x, 获取线搜上界
    '''

    t_min = float('inf')
    for i in range(p.shape[0]):
        if p[i] < 0:
            t_min = min(t_min, -x[i]/p[i])
    return t_min

    return np.max(A @ x - b)
def construct_p(B_k, N_k, r_N, I_B, I_ext, x):
    '''
    构造可行方向 p
    '''
    n = r_N.shape[0]
    p = np.zeros(B_k.shape[1]+N_k.shape[1])
    p_N = np.zeros(n)

    for i in range(n):
        if r_N[i] <= 0:
            p_N[i] = -r_N[i]
        else:
            p_N[i] = -x[I_ext[i]]*r_N[i]
    
    p_B = -np.linalg.inv(B_k) @ N_k @ p_N

    # 整合
    k = 0
    t = 0
    for i in range(B_k.shape[1]+N_k.shape[1]):
        if i in I_B:
            p[i] = p_B[k]
            k = k + 1
        else:
            p[i] = p_N[t]
            t = t + 1

    return p


def gradient(fun, x, epsilon=1e-10):
    '''
    求函数在 x 处的梯度
    '''
    dim = len(x)
    grad = np.zeros(dim)
    for i in range(dim):
        x_temp = x.copy()
        x_temp[i] += epsilon
        grad[i] = (fun(x_temp) - fun(x)) / epsilon
    return grad

def norm2(x):
    '''
    求向量的 L2 范数
    '''
    return np.sqrt(sum(x[i]**2 for i in range(len(x))))


def set_dynamic_fun(x, p):
    '''
    动态修改 fun 的参数
    '''
    def dynamic_fun(t):
        return Objective(x + t*p)
    return dynamic_fun

main()