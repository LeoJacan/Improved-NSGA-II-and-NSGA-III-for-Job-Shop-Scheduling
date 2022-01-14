# Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling
Learning how to implement a improved NSGA-II algorithm and NSGA-III algorithm for job shop scheduling problem in python.  
_2022/01/14_

本主題主要介紹如何透過基因演算法（Genetic Algorithm, GA）中的非凌越排序基因演算法（Nondominated Sorting Genetic Algorithm II, NSGA-II）之改良種演算法來求解 Job Shop 生產排程問題，以及另外使用第三代的非凌越排序基因演算法（NSGA-III）求解。之前已經有介紹 [GA](https://github.com/wurmen/Genetic-Algorithm-for-Job-Shop-Scheduling-and-NSGA-II/blob/master/introduction/GA/GA.md) 及 [NSGA-II](https://github.com/wurmen/Genetic-Algorithm-for-Job-Shop-Scheduling-and-NSGA-II/blob/master/introduction/NSGA-II/NSGA-II.md) 的專題，本專題是針對 GA 與 NSGA-II 改良部分以及 NSGA-III 進行說明，最後再透過 Python 來展示實作結果。  
我們改良 NSGA-II 的方法分為兩種，第一個是線性調整交配突變機制，第二個是線性調整選擇中的菁英機制。

[問題描述](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/%E5%95%8F%E9%A1%8C%E6%8F%8F%E8%BF%B0.md)  
GA 融合差分演算法：[文章](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/GA%20%E8%9E%8D%E5%90%88%E5%B7%AE%E5%88%86%E6%BC%94%E7%AE%97%E6%B3%95.md) / [code](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/NSGA%20code/%E5%B7%AE%E5%88%86%E6%BC%94%E7%AE%97%E6%B3%95%E5%8A%A0%E5%9F%BA%E5%9B%A0%E6%BC%94%E7%AE%97%E6%B3%95%E8%AA%BF%E6%95%B4.py) / import 模組 [NSGA-II](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/NSGA%20code/NSGAII.py)  
NSGA-II 線性調整　：[文章](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/NSGA-II%20%E7%B7%9A%E6%80%A7%E8%AA%BF%E6%95%B4.md) / [code1(交配突變)](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/NSGA%20code/NSGA-II%E4%BA%A4%E9%85%8D%E7%AA%81%E8%AE%8A%E7%B7%9A%E6%80%A7%E8%AA%BF%E6%95%B4.py) / [code2(菁英機制)](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/NSGA%20code/NSGA-II%E8%8F%81%E8%8B%B1%E7%B7%9A%E6%80%A7%E8%AA%BF%E6%95%B4.py) / import 模組 [NSGA-II](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/NSGA%20code/NSGAII.py)  
NSGA-III　　　　　：[文章](https://github.com/LeoJacan/Improved-NSGA-II-for-Job-Shop-Scheduling/blob/main/NSGA-III.md) / [code](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/NSGA%20code/untitled4.py) / import 模組 [NSGA-III](https://github.com/LeoJacan/Improved-NSGA-II-and-NSGA-III-for-Job-Shop-Scheduling/blob/main/NSGA%20code/NSGAIII.py)  
