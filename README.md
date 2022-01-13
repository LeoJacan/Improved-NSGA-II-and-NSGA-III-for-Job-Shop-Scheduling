# Improved-NSGA-II-for-Job-Shop-Scheduling
Learning how to implement a improved NSGA-II algorithm for job shop scheduling problem in python.

本主題主要介紹如何透過基因演算法（Genetic Algorithm, GA）中的非凌越排序基因演算法（Nondominated Sorting Genetic Algorithm II, NSGA-II）之改良種演算法來求解 Job Shop 排程問題。一開始會先進行 GA 及 NSGA-II 的概念介紹，再來針對 NSGA-II 改良部分進行說明，最後再透過 Python 來展示實作結果。

在半導體業中除了完成大量的訂單外，原料的節省不浪費也可以相對性降低不少成本。
利用多目標演算法平衡總完工時間（makespan）以及報廢量，Min-min是此實驗的目標將時間及報廢量最小化。
