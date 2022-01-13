### ✏前言

前面的文章有介紹了基因演算法 (GA) 衍伸出的非凌越排序基因演算法 (NSGA-II)，而本文將要介紹的 NSGA-III 是從 NSGA-II 的基礎進行改良改良而來，因此該演算法的架構與 NSGA-II 相似，NSGA-III 由 Kalyanmoy Deb 與 Himanshu Jain 於 2014 年發表，專門被用來求解高目標最佳化問題。本篇文章將要介紹何謂 NSGA-III ，並透過 Python 來進行實作，求解生產排程問題。

### ✏NSGA-III 的架構

 NSGA-II 與 NSGA-III 的差異之處就在於兩者為每個解(染色體)進行非凌越排序 (non-dominated sort) 後，NSGA-II 接下來是計算擁擠距離 (crowding-distance)，NSGA-III則是把此部分改成利基保存機制 (niche-preservation operation) ，流程架構可由下圖所示，紅框部分即是 NSGA-III 的利基保存機制，接下來也會針對此步驟進行說明。

![image](https://user-images.githubusercontent.com/97092223/149080169-d22cb9e2-9116-4813-a135-60cbfece4819.png)

### 🛠利基保存機制 (Niche-preservation operation)

 NSGA-III 利用利基保存機制來維持搜尋解的分散程度(diversity)，並加入一組參考點幫助每一代的染色體(解)演化。經過前一步的非凌越排序過後，從優先度最高的階層開始加入下一代的母代，最後一個有機會加入的階層利用利基保存機制挑選到下一代的解，為了幫助機制進行比較，藉由建構超平面(Hyperplane)來對解進行適應性正規化(adaptive normalization)，正規化後將選重的凌越階層中的所有解與參考方向進行關連，可由以下雙目標圖來舉例。
　　![image](https://user-images.githubusercontent.com/97092223/149094966-71017dc4-eaf8-4337-9609-b770fc19a0a4.png)
關聯後就到了最後的選擇機制，優先考慮每個參考方向上關連的解數量，判定利基計數的有或無來挑選，若一條參考方向上有離該參考方向垂直距離的唯一最小解，則優先挑選加入下一代母代，反之，若該參考方向沒有或超過一個最小垂直距離的解，則隨機挑選一個解到下一代。
利基保存機制的詳細計算步驟如下所示：

📎1.參考點的初始化：在目標數為 m 時，從 m-1 維的單形中均勻取樣產生參考點，得出參考點數 N，H 代表目標方向上的分割數  
![image](https://user-images.githubusercontent.com/97092223/149109622-f14fc019-d90a-49b8-9148-4cb88c0a27e3.png)

📎2.適應性解正規化：由解(S)裡面在各目標方向上最小的目標組成理想點(ideal point)(x)，接著從 S 中找出 m 個目標方向上的極限點(extreme points)，極限點為使 m 個目標方向的向量上形成標量化函數(achievement scalarizing function, ASF)值最小的解，ASF函數計算式如下：  
![image](https://user-images.githubusercontent.com/97092223/149116640-2cc26713-ef57-4130-aa6b-e9cddddbb72c.png)

📎3.產生超平面：由此 m 個極限店構成一個超平面，接著可以計算出此超平面在各目標方向上到原點的截距 a。有了以上資訊，便可將 S 裡的解 x 經由以下算式進行正規化：  
![image](https://user-images.githubusercontent.com/97092223/149118141-72b05010-97db-4b69-9a50-ec331ce73581.png)

📎4.關連：正規化之後，為了在環境選擇中讓參考點補助進行環境選擇，將選上之非凌越前緣階層中的解 S 與參考點 Z 進行關連，首先將參考點與原點連線形成參考線，接著計算每個解到最近參考線的垂直距離，若參考線擁有一個以上垂直距離最短的解，代表互相關連，同時利基計數便加一。

### 🛠選擇機制
經由非凌越前緣排序與上述過程，一代 Population 裡的所有染色體(解)經由非凌越前緣的排序階層(越小)優先挑選出來加入下一代親代，滿足親代數量的最後一階層前緣則是使用利基保存機制挑選剩餘的解，利基計數最小的參考方向優先判定，若該參考方向上有垂直距離最小的解則選擇之，若不只一個垂直最小解則隨機選擇加入下一代親代染色體，最終再次經過交配與突變產生新的 Population。
