# 代码 3 一维搜索方法
import numpy as np

def function(t):
    '''
    在此编辑函数形式
    min x(t)
    '''
    return t**3 - 2*t + 1

def golden_section_search(a, b, fun, epsilon=1e-5):
    '''
    0.618法（近似黄金分割法）
    '''

    while True:
        t1 = a + 0.382*(b-a)
        t2 = a + 0.618*(b-a)

        phi1 = fun(t1)
        phi2 = fun(t2)

        if phi1 <= phi2:
            if t2 - a <= epsilon:
                return t1
            else:
                b = t2
        else:
            if b - t1 <= epsilon:
                return t2
            else:
                a = t1

def newton_method(fun, x0, epsilon=1e-5):
    '''
    牛顿法
    '''
    x = x0

    while True:
        if abs(derivation(fun, x)) <= epsilon :
            return x
        elif derivation_2(fun, x) == 0:
            print("解体失败")
            return None
        else:
            x = x - derivation(fun, x)/derivation_2(fun, x)
            if abs(derivation(fun, x)/derivation_2(fun, x)) <= epsilon:
                return x

def derivation(fun, x, epsilon=1e-5):
    '''
    求导
    '''
    return (fun(x+epsilon)-fun(x))/epsilon

def derivation_2(fun, x, epsilon=1e-5):
    '''
    求二阶导
    '''
    return (fun(x+epsilon)-2*fun(x)+fun(x-epsilon))/epsilon**2

def get_search_bound(fun):
    '''
    划界算法
    确定初始搜索区间提高效率
    a: 下界
    b: 上界
    '''
    # 参数 
    h = 1
    lam = 2

    # 初始化
    a = 0
    b = h*lam
    c = 0

    # 循环
    while fun(c) > fun(b):
        a = c
        c = b
        h = lam*h
        b = b + h
    
    return a, b

# print("0.618法：")
# print(golden_section_search(0, 1, function))
# print("     ")
# print("牛顿法：")
# print(newton_method(function, 0.5))
# print("     ")