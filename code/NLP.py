# 代码 4：非线性规划——无约束
import numpy as np
import matplotlib.pyplot as plt
import LineSearch

def main():
    '''
    完成一些问题的设置
    和最后解的显示
    '''

    # 设置初始值
    x0 = np.array([3., 2.])

    # 设置精度
    epsilon = 1e-6

    # 调用最速下降法
    x_opt, x, y = steepest_descent_method(Objective, x0, epsilon)
    print("最速下降法：")
    print('x_opt =', x_opt)
    print('f(x_opt) =', Objective(x_opt))

    # 调用共轭梯度法
    x_opt1, x1, y1 = conjugate_direction_method(Objective, x0, epsilon)
    print("共轭方向法：")
    print('x_opt =', x_opt1)
    print('f(x_opt) =', Objective(x_opt1))

    # 最速下降图
    plt.figure()
    plt.plot(x, y, c='b', label='steepest descent', marker='*') 
    plt.legend(loc = 'lower center')
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',0))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',0))

    # 共轭方向图
    plt.figure()
    plt.plot(x1, y1, c='r', label='conjugate direction', marker='o')
    plt.legend(loc = 'lower center')
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',0))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',0))
    
    plt.show()

    
def Objective(x):
    return 1/3*x[0]**2 + 1/2*x[1]**2

def steepest_descent_method(fun, x0, epsilon=1e-5):
    '''
    最速下降法
    '''
    x = x0
    x_list = np.array(x[0])
    y_list = np.array(x[1])
    d = norm2(gradient(fun, x))
    k = 0
    while abs(d) > epsilon:
        
        p = -gradient(fun, x)

        t = 0
        dynamic_fun = set_dynamic_fun(x, p)
        a,b = LineSearch.functionget_search_bound(dynamic_fun)
        t = LineSearch.golden_section_search(a,b,dynamic_fun, epsilon)

        x = x + t*p
        k = k + 1

        x_list = np.append(x_list, x[0])
        y_list = np.append(y_list, x[1])

        print("第", k ,"次迭代")

        d = norm2(gradient(fun, x))

    return x, x_list, y_list

def conjugate_direction_method(fun, x0, epsilon=1e-5):
    '''
    共轭方向法
    '''
    # step 1
    x = x0
    x_list = np.array(x[0])
    y_list = np.array(x[1])

    # step 2
    d = norm2(gradient(fun, x))
    k = 0

    # step 3
    p = -gradient(fun, x)

    while d > epsilon:
        print("第", k+1 ,"次迭代")

        # step 4
        dynamic_fun = set_dynamic_fun(x, p)
        a,b = LineSearch.get_search_bound(dynamic_fun)
        t = LineSearch.golden_section_search(a,b,dynamic_fun, epsilon)
        x_k = x + t*p
        x_list = np.append(x_list, x_k[0])
        y_list = np.append(y_list, x_k[1])

        # step 5
        d = norm2(gradient(fun, x_k))
        if d <= epsilon:
            return x_k, x_list, y_list

        # step 6
        if k+1 == x.shape[0]:
            continue
        else:
            # step 7
            # F-R 公式 
            lam_k = norm2(gradient(fun,x_k))**2/norm2(gradient(fun, x))**2
            p = -gradient(fun, x_k) + lam_k*p

            k = k + 1
            x = x_k



    return x, x_list, y_list
 

    

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
    def dynamic_fun(t):
        return Objective(x + t*p)
    return dynamic_fun


main()