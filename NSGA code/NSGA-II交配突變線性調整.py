# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 15:14:28 2022

@author: lian
"""


from pandas import Series
import pandas as pd
import numpy as np
import random
import copy
import time

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# import Draw_GanttChart as Draw
from datetime import datetime, timedelta
import plotly.express as px
from collections import Counter
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
from plotly.figure_factory import create_gantt
from operator import itemgetter
import datetime
import plotly as py
import time
import xlrd
All_Machines_Name = ['1','2','3','4','5','6','7','8','9','10']
#跑資料=========================================================================
test_arry=[]
Num_job=100
Num_Machine=10
Num_chromosome=100
a = pd.read_excel("semiconductor_data.xlsx",sheet_name="WIP")
b = pd.read_excel("semiconductor_data.xlsx",sheet_name="EQP_RECIPE")
c=pd.read_excel("semiconductor_data.xlsx",sheet_name="set_up_time")
WIP_data = copy.deepcopy(a)
EQP_RECIPE_data = copy.deepcopy(b)
setup_time=copy.deepcopy(c)
initialization_chromosome=np.zeros((Num_chromosome*2,Num_job*2))
select_machine=np.zeros((Num_job,11))
cross_rate=0.8#交配率
mutation_rate=0.2#突變率
lenght=np.arange(Num_chromosome)
iteration=500
iteration_makespan=[]
iteration_Scrap=[]


cross_rate_Variety=np.linspace(0.5,1,num = 500,endpoint=True)
RATE=np.linspace(0.2,1,num = 500,endpoint=True)
for i in range(1,11):
    if i ==10:
        CAN=WIP_data["CANRUN_TOOL"].str.split('EQP0',10,expand=True)
        CANRUN_TOOL_data = copy.deepcopy(CAN)
    else:
        CAN=WIP_data["CANRUN_TOOL"].str.split('EQP0',i,expand=True)
        CANRUN_TOOL_data = copy.deepcopy(CAN)
CANRUN_TOOL_data=CANRUN_TOOL_data.fillna(0)

#編碼#最初染色體
for k in range(1,iteration+1):
    print("第 %d 帶" % k)
    start_time =time.time()
    if k == 1:
        for i in range(Num_chromosome):
            for j in range(Num_job*2):
                initialization_chromosome[i][j]=np.random.rand()
            
        #交配===========================================================================
        random.shuffle(lenght)
        parents_A =lenght[0:int(len(lenght)/2)]#分成A組
        parents_B =lenght[int(len(lenght)/2):]#分成B組
        for i in range(int((Num_chromosome*cross_rate_Variety[k-1])/2)):
            j=np.array((random.sample(range(Num_job*2),2)))#隨機抽選要交配之行數
            sorted(j)#由小排到大排區間
            children_A=initialization_chromosome[parents_A[i]][j[0]:j[1]]#尋找要更換的染色體位置
            children_B=initialization_chromosome[parents_B[i]][j[0]:j[1]]
            initialization_chromosome[Num_chromosome+2*i,0:Num_job*2]=initialization_chromosome[parents_A[i]][0:Num_job*2]#先將要交配之染色體從第10條開始往下存取
            initialization_chromosome[(Num_chromosome+1)+2*i,0:Num_job*2]=initialization_chromosome[parents_B[i]][0:Num_job*2]
            initialization_chromosome[Num_chromosome+2*i,j[0]:j[1]]=children_B#進行交配
            initialization_chromosome[(Num_chromosome+1)+2*i,j[0]:j[1]]=children_A
                            #交配完成
        #突變============================================================================
            mutation_rows=0
            lenght2=np.arange(Num_chromosome+(Num_chromosome*cross_rate_Variety[k-1]))
            random.shuffle(lenght2)
            i+=1
            parents_C = lenght2[0:int((Num_chromosome*(1-cross_rate_Variety[k-1]))/2)]
            parents_D = lenght2[int((Num_chromosome*(1-cross_rate_Variety[k-1]))/2):]
            for x in range(0,int((Num_chromosome*(1-cross_rate_Variety[k-1]))/2)):
                initialization_chromosome[Num_chromosome+2*i,0:Num_job*2]=initialization_chromosome[int(parents_C[x])][0:Num_job*2]
                initialization_chromosome[(Num_chromosome+1)+2*i,0:Num_job*2]=initialization_chromosome[int(parents_D[x])][0:Num_job*2]
                mutation_rows1=np.array((random.sample(range(Num_job*2),1)))#隨機抽選要突變之行數1
                mutation_rows1_1=np.array((random.sample(range(Num_job*2),1)))#隨機抽選要突變之行數1-1
                mutation_randompoint=np.random.rand()#單點突變數1
                mutation_randompoint1_1=np.random.rand()#單點突變數1-1
                initialization_chromosome[Num_chromosome+2*i][mutation_rows1]=mutation_randompoint
                initialization_chromosome[Num_chromosome+2*i][mutation_rows1_1]=mutation_randompoint1_1
                mutation_rows2=np.array((random.sample(range(Num_job*2),1)))#隨機抽選要突變之行數2
                mutation_rows2_1=np.array((random.sample(range(Num_job*2),1)))#隨機抽選要突變之行數2-1
                mutation_randompoint2=np.random.rand()#單點突變數2
                mutation_randompoint2_1=np.random.rand()#單點突變數2-1
                initialization_chromosome[(Num_chromosome+1)+2*i][mutation_rows2]=mutation_randompoint2
                initialization_chromosome[(Num_chromosome+1)+2*i][mutation_rows2_1]=mutation_randompoint2_1
                i+=1
        #選機============================================================================
        count_canuse_machine=0
        canuse_machine=[]
        for i in range(Num_job):
            for  j in range(1,Num_Machine):
                if int(CANRUN_TOOL_data.iloc[i][j])>0:
                    count_canuse_machine=count_canuse_machine+1
            canuse_machine.append(count_canuse_machine)
            count_canuse_machine=0
        weight_number=1
        All_select_machine=[]  
        for z in range(Num_chromosome*2):
            weight_number=1
            for i in range(Num_job):
                for j in range(1,Num_Machine):
                    division_number=1/canuse_machine[i]
                    if division_number*weight_number<initialization_chromosome[z][i]:
                        weight_number+=1
                    elif division_number*weight_number>initialization_chromosome[z][i]:
                        break;
                #print(i,weight_number,CANRUN_TOOL_data.iloc[i][weight_number])
                select_machine[i][0]=i
                select_machine[i][1]=weight_number
                select_machine[i][2]=int(CANRUN_TOOL_data.iloc[i][weight_number])
                select_machine[i][3]=initialization_chromosome[z][100+i]
                weight_number=1
            select_machine=pd.DataFrame(select_machine)
            #將蒐集成的選機及排序加上標籤及做成dataframe
            select_machine2={'工作':Series(select_machine.iloc[0:100,0]),
                                 'b':Series(select_machine.iloc[0:100,1]),
                                 '機台':Series(select_machine.iloc[0:100,2]),
                                 '選機隨機亂數':Series(select_machine.iloc[0:100,3]),
                                 '運作時間':Series(select_machine.iloc[0:100,4]),
                                 '起始時間':Series(select_machine.iloc[0:100,5]),
                                 '終止時間':Series(select_machine.iloc[0:100,6]),
                                 '加換模起始時間':Series(select_machine.iloc[0:100,7]),
                                 '加換模終止時間':Series(select_machine.iloc[0:100,8]),
                                 '換模及延誤起始時間':Series(select_machine.iloc[0:100,9]),
                                 '換模及延誤終止時間':Series(select_machine.iloc[0:100,10]),
                                 }
            select_machine_finish=pd.DataFrame(select_machine2)
            #print(select_machine_finish)   
            select_machine_finish.sort_values(by=['機台' ,'選機隨機亂數'], ascending=(True,False) ,inplace=True)
            All_select_machine.append(select_machine_finish)
            select_machine=select_machine.values#最後再把dataframe轉回去給下一個染色體使用
        
        for z in range(Num_chromosome*2):
            for  i in range(100):
                #X=All_select_machine[z].iloc[i][2]==EQP_RECIPE_data["EQP_ID"]
                start=int(All_select_machine[z].iloc[i][2])
                for j in range((start-1)*10,(start-1)*10+10):
                        if WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][8]==EQP_RECIPE_data.iloc[j][1]:
                            #print(int(EQP_RECIPE_data.iloc[j][2]))
                            All_select_machine[z].iloc[i][4]=int(EQP_RECIPE_data.iloc[j][2])
        #片數時間計算
        for z in range(Num_chromosome*2):
            for  i in range(0,100):
                All_select_machine[z].iloc[i][4]=All_select_machine[z].iloc[i][4]*(WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][3]/25)
        #計算時間
        for z in range(Num_chromosome*2):
            for  i in range(1,100):
                if All_select_machine[z].iloc[i][2]!=All_select_machine[z].iloc[i-1][2]:
                    All_select_machine[z].iloc[i][6]=All_select_machine[z].iloc[i][5]+All_select_machine[z].iloc[i][4]
        
                    All_select_machine[z].iloc[i][5]=0
                    #All_select_machine[z].iloc[i][7]=0##
                elif All_select_machine[z].iloc[i][2]==All_select_machine[z].iloc[i-1][2]:
                    All_select_machine[z].iloc[0][6]=All_select_machine[z].iloc[0][5]+All_select_machine[z].iloc[0][4]
                    All_select_machine[z].iloc[i][6]=All_select_machine[z].iloc[i-1][6]+All_select_machine[z].iloc[i][4]
                    All_select_machine[z].iloc[i][5]=All_select_machine[z].iloc[i-1][6]
        
        #整備時間            
        for z in range(Num_chromosome*2):
            for  i in range(0,99):
                if All_select_machine[z].iloc[i+1][2]!=All_select_machine[z].iloc[i][2]:
                    All_select_machine[z].iloc[i+1][7]=0##
                    All_select_machine[z].iloc[i+1][8]=All_select_machine[z].iloc[i+1][7]+All_select_machine[z].iloc[i+1][4]##
                    
                elif All_select_machine[z].iloc[i+1][2]==All_select_machine[z].iloc[i][2]:
                    setuptime_data=setup_time.iloc[int(All_select_machine[z].iloc[i][0])][int(All_select_machine[z].iloc[i+1][0])+1]##
                    All_select_machine[z].iloc[0][8]=All_select_machine[z].iloc[0][4]+All_select_machine[z].iloc[0][5]
                    All_select_machine[z].iloc[i+1][7]=All_select_machine[z].iloc[i][8]+int(setuptime_data)##
                    All_select_machine[z].iloc[i+1][8]=All_select_machine[z].iloc[i+1][7]+All_select_machine[z].iloc[i+1][4]##
        #延誤時間
        for z in range(Num_chromosome*2):
            if WIP_data.iloc[int(All_select_machine[z].iloc[0][0])][5] > All_select_machine[0].iloc[0][7]:
                All_select_machine[z].iloc[0][9]=WIP_data.iloc[int(All_select_machine[z].iloc[0][0])][5]
                All_select_machine[z].iloc[0][10]=All_select_machine[z].iloc[0][9]+All_select_machine[z].iloc[0][4]
            else:
                All_select_machine[z].iloc[0][9]=All_select_machine[z].iloc[0][7]
                All_select_machine[z].iloc[0][10]=All_select_machine[z].iloc[0][8]
            for i in range(1,100):
                gap = All_select_machine[z].iloc[i][8]-All_select_machine[z].iloc[i][7]
                Next_generation_gap=All_select_machine[z].iloc[i][7]-All_select_machine[z].iloc[i-1][8]
                if All_select_machine[z].iloc[i][7]==0 and WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5]>All_select_machine[z].iloc[i][7]:
                    All_select_machine[z].iloc[i][9]=WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5]
                    All_select_machine[z].iloc[i][10]=All_select_machine[z].iloc[i][9]+gap
                if All_select_machine[z].iloc[i][7]==0 and WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5]<=All_select_machine[z].iloc[i][7]:
                    All_select_machine[z].iloc[i][9]=All_select_machine[z].iloc[i][7]
                    All_select_machine[z].iloc[i][10]=All_select_machine[z].iloc[i][9]+gap
                
                if WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5] > All_select_machine[z].iloc[i-1][10]+Next_generation_gap and  All_select_machine[z].iloc[i][7]!=0:
                    All_select_machine[z].iloc[i][9]=WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5]
                    All_select_machine[z].iloc[i][10]=All_select_machine[z].iloc[i][9]+gap
                if WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5] <= All_select_machine[z].iloc[i-1][10]+Next_generation_gap and  All_select_machine[z].iloc[i][7]!=0:
                    All_select_machine[z].iloc[i][9]=All_select_machine[z].iloc[i-1][10]+Next_generation_gap
                    All_select_machine[z].iloc[i][10]=All_select_machine[z].iloc[i][9]+gap
                
        collect_chromosome_mqakespan_saveplace=np.zeros((Num_chromosome*2,4))
        scrapped_object_time_all=[]
        scrapped_object_time=np.zeros((Num_job,3))
        for z in range(Num_chromosome*2):
            for i in range(Num_job):
                 scrapped_object_time[i][0]=All_select_machine[z].iloc[i][9]
                 scrapped_object_time[i][1]=WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][4]*60
                 if scrapped_object_time[i][1]>scrapped_object_time[i][0]:
                     scrapped_object_time[i][2]=0
                 elif scrapped_object_time[i][1]<scrapped_object_time[i][0]:
                     scrapped_object_time[i][2]=1
            scrapped_object_time_all.append(scrapped_object_time)
            collect_chromosome_mqakespan_saveplace[z][2]=sum(scrapped_object_time_all[z][0:,2])
            scrapped_object_time=np.zeros((Num_job,3))
                    
        for z in range(Num_chromosome*2):
            a=max(All_select_machine[z].iloc[0:,10])
            collect_chromosome_mqakespan_saveplace[z][0]=a
            collect_chromosome_mqakespan_saveplace[z][1]=z
            collect_chromosome_mqakespan_saveplace[z][3]=collect_chromosome_mqakespan_saveplace[z][0]+collect_chromosome_mqakespan_saveplace[z][2]
        collect_chromosome_mqakespan_saveplace=pd.DataFrame(collect_chromosome_mqakespan_saveplace)
        collect_chromosome_mqakespan_s={'完工時間':Series(collect_chromosome_mqakespan_saveplace.iloc[0:Num_chromosome*2,0]),
                                      '工作':Series(collect_chromosome_mqakespan_saveplace.iloc[0:Num_chromosome*2,1]),
                                      '報廢量':Series(collect_chromosome_mqakespan_saveplace.iloc[0:Num_chromosome*2,2]),
                                      '目標值':Series(collect_chromosome_mqakespan_saveplace.iloc[0:Num_chromosome*2,3])
                                      }
        collect_chromosome_mqakespan=pd.DataFrame(collect_chromosome_mqakespan_s)
        collect_chromosome_mqakespan.sort_values(by=['完工時間'], ascending=(True) ,inplace=True)
        
        # NSGA_TEST = collect_chromosome_mqakespan.iloc[:,[0,2]]
        # import pymoo
        # import NSGAII
        # NSGAII.ValuesCal(NSGA_TEST,200)
        # test_NSGA_TEST = NSGAII.RankforNSGAII(NSGA_TEST)
        # print(test_NSGA_TEST.iloc[0])

        
        
        
        #挑選
        Elite_rate=RATE[k-1]#菁英制
        Level_Roulette=0.8#等級輪盤
        

        iteration_makespan.append(collect_chromosome_mqakespan.iloc[0][0])
        
        NSGA_TEST = collect_chromosome_mqakespan.iloc[:,[0,2]]
        import pymoo
        import NSGAII
        NSGAII.ValuesCal(NSGA_TEST,200)
        test_NSGA_TEST = NSGAII.RankforNSGAII(NSGA_TEST)
        print(test_NSGA_TEST.iloc[0])
        initialization_chromosome_select_nb=[]
        initialization_chromosome_select=np.zeros((Num_chromosome,Num_job*2))#新的母代
        test_NSGA_TEST_2=test_NSGA_TEST.sort_values("Rank", ascending=True)
        for i in range(0,int(Num_chromosome*Elite_rate)):
                initialization_chromosome_select_nb.append(test_NSGA_TEST_2.index[i])
                initialization_chromosome_select[i][0:200]=initialization_chromosome[test_NSGA_TEST.index[i],0:200]
        drowlot_nb=np.random.rand(100-int(Num_chromosome*Elite_rate))
        #sumnb=sum(test_NSGA_TEST_2.iloc[int(Num_chromosome*Elite_rate):,4])
        nsga_ranknb=np.zeros((Num_chromosome*2-int(Num_chromosome*Elite_rate),4))
        x=1
        for i in range(Num_chromosome*2-1,int(Num_chromosome*Elite_rate)-1,-1):
            if test_NSGA_TEST_2.iloc[i][4]!=test_NSGA_TEST_2.iloc[i-1][4]:
                nsga_ranknb[i-int(Num_chromosome*Elite_rate),0]=x
                x+=1
            elif test_NSGA_TEST_2.iloc[i][4]==test_NSGA_TEST_2.iloc[i-1][4]:
                nsga_ranknb[i-int(Num_chromosome*Elite_rate),0]=x
        sumnb=sum(nsga_ranknb[:,0])
        for i in range(0,int(Num_chromosome*2)-int(Num_chromosome*Elite_rate)):
            nsga_ranknb[i][1]=nsga_ranknb[i][0]/sumnb
            nsga_ranknb[0][2]=nsga_ranknb[0][1]
            nsga_ranknb[i][3]=test_NSGA_TEST_2.index[i+int(Num_chromosome*Elite_rate)]
        for i in range(1,int(Num_chromosome*2)-int(Num_chromosome*Elite_rate)):
            nsga_ranknb[i][2]=nsga_ranknb[i-1][2]+nsga_ranknb[i][1]
        for i in range(len(drowlot_nb)):
            for j in range(len(nsga_ranknb[:,2])):
                if drowlot_nb[i]<nsga_ranknb[j][2]:
                    print("選第",nsga_ranknb[j][3],"個")
                    initialization_chromosome_select[i+int(Num_chromosome*Elite_rate)][0:200]=initialization_chromosome[[int(nsga_ranknb[j][3])],0:200]
                    break
        iteration_Scrap.append(test_NSGA_TEST.iloc[0][1])

 
        # iteration_makespan.append(collect_chromosome_mqakespan.iloc[0][0])
    elif k > 1:
        initialization_chromosome=np.zeros((Num_chromosome*2,Num_job*2))
        for i in range(0,Num_chromosome):
            for j in range(0,Num_job*2):
                initialization_chromosome[i][j]=initialization_chromosome_select[i][j]                
        random.shuffle(lenght)
        parents_A =lenght[0:int(len(lenght)/2)]#分成A組
        parents_B =lenght[int(len(lenght)/2):]#分成B組
        for i in range(int((Num_chromosome*cross_rate_Variety[k-1])/2)):
            j=np.array((random.sample(range(Num_job*2),2)))#隨機抽選要交配之行數
            sorted(j)#由小排到大排區間
            children_A=initialization_chromosome[parents_A[i]][j[0]:j[1]]#尋找要更換的染色體位置
            children_B=initialization_chromosome[parents_B[i]][j[0]:j[1]]
            initialization_chromosome[Num_chromosome+2*i,0:Num_job*2]=initialization_chromosome[parents_A[i]][0:Num_job*2]#先將要交配之染色體從第10條開始往下存取
            initialization_chromosome[(Num_chromosome+1)+2*i,0:Num_job*2]=initialization_chromosome[parents_B[i]][0:Num_job*2]
            initialization_chromosome[Num_chromosome+2*i,j[0]:j[1]]=children_B#進行交配
            initialization_chromosome[(Num_chromosome+1)+2*i,j[0]:j[1]]=children_A
                            #交配完成
        #突變============================================================================
            mutation_rows=0
            lenght2=np.arange(Num_chromosome+(Num_chromosome*cross_rate_Variety[k-1]))
            random.shuffle(lenght2)
            i+=1
            parents_C = lenght2[0:int((Num_chromosome*(1-cross_rate_Variety[k-1]))/2)]
            parents_D = lenght2[int((Num_chromosome*(1-cross_rate_Variety[k-1]))/2):]
            for x in range(0,int((Num_chromosome*(1-cross_rate_Variety[k-1]))/2)):
                initialization_chromosome[Num_chromosome+2*i,0:Num_job*2]=initialization_chromosome[int(parents_C[x])][0:Num_job*2]
                initialization_chromosome[(Num_chromosome+1)+2*i,0:Num_job*2]=initialization_chromosome[int(parents_D[x])][0:Num_job*2]
                mutation_rows1=np.array((random.sample(range(Num_job*2),1)))#隨機抽選要突變之行數1
                mutation_rows1_1=np.array((random.sample(range(Num_job*2),1)))#隨機抽選要突變之行數1-1
                mutation_randompoint=np.random.rand()#單點突變數1
                mutation_randompoint1_1=np.random.rand()#單點突變數1-1
                initialization_chromosome[Num_chromosome+2*i][mutation_rows1]=mutation_randompoint
                initialization_chromosome[Num_chromosome+2*i][mutation_rows1_1]=mutation_randompoint1_1
                mutation_rows2=np.array((random.sample(range(Num_job*2),1)))#隨機抽選要突變之行數2
                mutation_rows2_1=np.array((random.sample(range(Num_job*2),1)))#隨機抽選要突變之行數2-1
                mutation_randompoint2=np.random.rand()#單點突變數2
                mutation_randompoint2_1=np.random.rand()#單點突變數2-1
                initialization_chromosome[(Num_chromosome+1)+2*i][mutation_rows2]=mutation_randompoint2
                initialization_chromosome[(Num_chromosome+1)+2*i][mutation_rows2_1]=mutation_randompoint2_1
                i+=1
        #選機============================================================================
        count_canuse_machine=0
        canuse_machine=[]
        for i in range(Num_job):
            for  j in range(1,Num_Machine):
                if int(CANRUN_TOOL_data.iloc[i][j])>0:
                    count_canuse_machine=count_canuse_machine+1
            canuse_machine.append(count_canuse_machine)
            count_canuse_machine=0
        weight_number=1
        All_select_machine=[]  
        for z in range(Num_chromosome*2):
            weight_number=1
            for i in range(Num_job):
                for j in range(1,Num_Machine):
                    division_number=1/canuse_machine[i]
                    if division_number*weight_number<initialization_chromosome[z][i]:
                        weight_number+=1
                    elif division_number*weight_number>initialization_chromosome[z][i]:
                        break;
                #print(i,weight_number,CANRUN_TOOL_data.iloc[i][weight_number])
                select_machine[i][0]=i
                select_machine[i][1]=weight_number
                select_machine[i][2]=int(CANRUN_TOOL_data.iloc[i][weight_number])
                select_machine[i][3]=initialization_chromosome[z][100+i]
                weight_number=1
            select_machine=pd.DataFrame(select_machine)
            #將蒐集成的選機及排序加上標籤及做成dataframe
            select_machine2={'工作':Series(select_machine.iloc[0:100,0]),
                                 'b':Series(select_machine.iloc[0:100,1]),
                                 '機台':Series(select_machine.iloc[0:100,2]),
                                 '選機隨機亂數':Series(select_machine.iloc[0:100,3]),
                                 '運作時間':Series(select_machine.iloc[0:100,4]),
                                 '起始時間':Series(select_machine.iloc[0:100,5]),
                                 '終止時間':Series(select_machine.iloc[0:100,6]),
                                 '加換模起始時間':Series(select_machine.iloc[0:100,7]),
                                 '加換模終止時間':Series(select_machine.iloc[0:100,8]),
                                 '換模及延誤起始時間':Series(select_machine.iloc[0:100,9]),
                                 '換模及延誤終止時間':Series(select_machine.iloc[0:100,10]),
                                 }
            select_machine_finish=pd.DataFrame(select_machine2)
            #print(select_machine_finish)   
            select_machine_finish.sort_values(by=['機台' ,'選機隨機亂數'], ascending=(True,False) ,inplace=True)
            All_select_machine.append(select_machine_finish)
            select_machine=select_machine.values#最後再把dataframe轉回去給下一個染色體使用
        
        for z in range(Num_chromosome*2):
            for  i in range(100):
                #X=All_select_machine[z].iloc[i][2]==EQP_RECIPE_data["EQP_ID"]
                start=int(All_select_machine[z].iloc[i][2])
                for j in range((start-1)*10,(start-1)*10+10):
                        if WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][8]==EQP_RECIPE_data.iloc[j][1]:
                            #print(int(EQP_RECIPE_data.iloc[j][2]))
                            All_select_machine[z].iloc[i][4]=int(EQP_RECIPE_data.iloc[j][2])
        #片數時間計算
        for z in range(Num_chromosome*2):
            for  i in range(0,100):
                All_select_machine[z].iloc[i][4]=All_select_machine[z].iloc[i][4]*(WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][3]/25)
        #計算時間
        for z in range(Num_chromosome*2):
            for  i in range(1,100):
                if All_select_machine[z].iloc[i][2]!=All_select_machine[z].iloc[i-1][2]:
                    All_select_machine[z].iloc[i][6]=All_select_machine[z].iloc[i][5]+All_select_machine[z].iloc[i][4]
        
                    All_select_machine[z].iloc[i][5]=0
                    #All_select_machine[z].iloc[i][7]=0##
                elif All_select_machine[z].iloc[i][2]==All_select_machine[z].iloc[i-1][2]:
                    All_select_machine[z].iloc[0][6]=All_select_machine[z].iloc[0][5]+All_select_machine[z].iloc[0][4]
                    All_select_machine[z].iloc[i][6]=All_select_machine[z].iloc[i-1][6]+All_select_machine[z].iloc[i][4]
                    All_select_machine[z].iloc[i][5]=All_select_machine[z].iloc[i-1][6]
        
        #整備時間            
        for z in range(Num_chromosome*2):
            for  i in range(0,99):
                if All_select_machine[z].iloc[i+1][2]!=All_select_machine[z].iloc[i][2]:
                    All_select_machine[z].iloc[i+1][7]=0##
                    All_select_machine[z].iloc[i+1][8]=All_select_machine[z].iloc[i+1][7]+All_select_machine[z].iloc[i+1][4]##
                    
                elif All_select_machine[z].iloc[i+1][2]==All_select_machine[z].iloc[i][2]:
                    setuptime_data=setup_time.iloc[int(All_select_machine[z].iloc[i][0])][int(All_select_machine[z].iloc[i+1][0])+1]##
                    All_select_machine[z].iloc[0][8]=All_select_machine[z].iloc[0][4]+All_select_machine[z].iloc[0][5]
                    All_select_machine[z].iloc[i+1][7]=All_select_machine[z].iloc[i][8]+int(setuptime_data)##
                    All_select_machine[z].iloc[i+1][8]=All_select_machine[z].iloc[i+1][7]+All_select_machine[z].iloc[i+1][4]##
        #延誤時間
        for z in range(Num_chromosome*2):
            if WIP_data.iloc[int(All_select_machine[z].iloc[0][0])][5] > All_select_machine[0].iloc[0][7]:
                All_select_machine[z].iloc[0][9]=WIP_data.iloc[int(All_select_machine[z].iloc[0][0])][5]
                All_select_machine[z].iloc[0][10]=All_select_machine[z].iloc[0][9]+All_select_machine[z].iloc[0][4]
            else:
                All_select_machine[z].iloc[0][9]=All_select_machine[z].iloc[0][7]
                All_select_machine[z].iloc[0][10]=All_select_machine[z].iloc[0][8]
            for i in range(1,100):
                gap = All_select_machine[z].iloc[i][8]-All_select_machine[z].iloc[i][7]
                Next_generation_gap=All_select_machine[z].iloc[i][7]-All_select_machine[z].iloc[i-1][8]
                if WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5] > All_select_machine[z].iloc[i-1][10]+Next_generation_gap and  All_select_machine[z].iloc[i][7]!=0:
                    All_select_machine[z].iloc[i][9]=WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5]
                    All_select_machine[z].iloc[i][10]=All_select_machine[z].iloc[i][9]+gap
                elif WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5] <= All_select_machine[z].iloc[i-1][10]+Next_generation_gap and  All_select_machine[z].iloc[i][7]!=0:
                    All_select_machine[z].iloc[i][9]=All_select_machine[z].iloc[i-1][10]+Next_generation_gap
                    All_select_machine[z].iloc[i][10]=All_select_machine[z].iloc[i][9]+gap
                elif All_select_machine[z].iloc[i][7]==0 and WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5]>All_select_machine[z].iloc[i][7]:
                    All_select_machine[z].iloc[i][9]=WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5]
                    All_select_machine[z].iloc[i][10]=All_select_machine[z].iloc[i][9]+gap
                elif All_select_machine[z].iloc[i][7]==0 and WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][5]<=All_select_machine[z].iloc[i][7]:
                    All_select_machine[z].iloc[i][9]=All_select_machine[z].iloc[i][7]
                    All_select_machine[z].iloc[i][10]=All_select_machine[z].iloc[i][9]+gap
        collect_chromosome_mqakespan_saveplace=np.zeros((Num_chromosome*2,4))
        scrapped_object_time_all=[]
        scrapped_object_time=np.zeros((Num_job,3))
        for z in range(Num_chromosome*2):
            for i in range(Num_job):
                 scrapped_object_time[i][0]=All_select_machine[z].iloc[i][9]
                 scrapped_object_time[i][1]=WIP_data.iloc[int(All_select_machine[z].iloc[i][0])][4]*60
                 if scrapped_object_time[i][1]>scrapped_object_time[i][0]:
                     scrapped_object_time[i][2]=0
                 elif scrapped_object_time[i][1]<scrapped_object_time[i][0]:
                     scrapped_object_time[i][2]=1
            scrapped_object_time_all.append(scrapped_object_time)
            collect_chromosome_mqakespan_saveplace[z][2]=sum(scrapped_object_time_all[z][0:,2])
            scrapped_object_time=np.zeros((Num_job,3))
                    
        for z in range(Num_chromosome*2):
            a=max(All_select_machine[z].iloc[0:,10])
            collect_chromosome_mqakespan_saveplace[z][0]=a
            collect_chromosome_mqakespan_saveplace[z][1]=z
            collect_chromosome_mqakespan_saveplace[z][3]=collect_chromosome_mqakespan_saveplace[z][0]+collect_chromosome_mqakespan_saveplace[z][2]
        collect_chromosome_mqakespan_saveplace=pd.DataFrame(collect_chromosome_mqakespan_saveplace)
        collect_chromosome_mqakespan_s={'完工時間':Series(collect_chromosome_mqakespan_saveplace.iloc[0:Num_chromosome*2,0]),
                                      '工作':Series(collect_chromosome_mqakespan_saveplace.iloc[0:Num_chromosome*2,1]),
                                      '報廢量':Series(collect_chromosome_mqakespan_saveplace.iloc[0:Num_chromosome*2,2]),
                                      '目標值':Series(collect_chromosome_mqakespan_saveplace.iloc[0:Num_chromosome*2,3])
                                      }
        collect_chromosome_mqakespan=pd.DataFrame(collect_chromosome_mqakespan_s)
        collect_chromosome_mqakespan.sort_values(by=['完工時間'], ascending=(True) ,inplace=True)
        
        #挑選
        Elite_rate=RATE[k-1]#菁英制
        Level_Roulette=0.8#等級輪盤
        

        iteration_makespan.append(collect_chromosome_mqakespan.iloc[0][0])
        NSGA_TEST = collect_chromosome_mqakespan.iloc[:,[0,2]]
        import pymoo
        import NSGAII
        NSGAII.ValuesCal(NSGA_TEST,200)
        test_NSGA_TEST = NSGAII.RankforNSGAII(NSGA_TEST)
        print(test_NSGA_TEST.iloc[0])
        initialization_chromosome_select_nb=[]
        initialization_chromosome_select=np.zeros((Num_chromosome,Num_job*2))#新的母代
        test_NSGA_TEST_2=test_NSGA_TEST.sort_values("Rank", ascending=True)
        for i in range(0,int(Num_chromosome*Elite_rate)):
                initialization_chromosome_select_nb.append(test_NSGA_TEST_2.index[i])
                initialization_chromosome_select[i][0:200]=initialization_chromosome[test_NSGA_TEST.index[i],0:200]
        drowlot_nb=np.random.rand(100-int(Num_chromosome*Elite_rate))
        #sumnb=sum(test_NSGA_TEST_2.iloc[int(Num_chromosome*Elite_rate):,4])
        nsga_ranknb=np.zeros((Num_chromosome*2-int(Num_chromosome*Elite_rate),4))
        x=1
        for i in range(Num_chromosome*2-1,int(Num_chromosome*Elite_rate)-1,-1):
            if test_NSGA_TEST_2.iloc[i][4]!=test_NSGA_TEST_2.iloc[i-1][4]:
                nsga_ranknb[i-int(Num_chromosome*Elite_rate),0]=x
                x+=1
            elif test_NSGA_TEST_2.iloc[i][4]==test_NSGA_TEST_2.iloc[i-1][4]:
                nsga_ranknb[i-int(Num_chromosome*Elite_rate),0]=x
        sumnb=sum(nsga_ranknb[:,0])
        for i in range(0,int(Num_chromosome*2)-int(Num_chromosome*Elite_rate)):
            nsga_ranknb[i][1]=nsga_ranknb[i][0]/sumnb
            nsga_ranknb[0][2]=nsga_ranknb[0][1]
            nsga_ranknb[i][3]=test_NSGA_TEST_2.index[i+int(Num_chromosome*Elite_rate)]
        for i in range(1,int(Num_chromosome*2)-int(Num_chromosome*Elite_rate)):
            nsga_ranknb[i][2]=nsga_ranknb[i-1][2]+nsga_ranknb[i][1]
        for i in range(len(drowlot_nb)):
            for j in range(len(nsga_ranknb[:,2])):
                if drowlot_nb[i]<nsga_ranknb[j][2]:
                    print("選第",nsga_ranknb[j][3],"個")
                    initialization_chromosome_select[i+int(Num_chromosome*Elite_rate)][0:200]=initialization_chromosome[[int(nsga_ranknb[j][3])],0:200]
                    break
        iteration_Scrap.append(test_NSGA_TEST.iloc[0][1])

    end=time.time()
    print(end-start_time)
Gantt_chart_number=np.zeros((Num_job,4))
for  i in range(Num_job):
    Gantt_chart_number[i][0]=All_select_machine[int(collect_chromosome_mqakespan.iloc[0][1])].iloc[i][0]
    Gantt_chart_number[i][1]=All_select_machine[int(collect_chromosome_mqakespan.iloc[0][1])].iloc[i][9]
    Gantt_chart_number[i][2]=All_select_machine[int(collect_chromosome_mqakespan.iloc[0][1])].iloc[i][10]
    Gantt_chart_number[i][3]=All_select_machine[int(collect_chromosome_mqakespan.iloc[0][1])].iloc[i][2]
    
df = []
def draw_machine(index,value):
    '''把Gantt_chart_number 結束時間排序過後的第一條拿來畫圖!!'''
    for i in range(len(Gantt_chart_number)):
        if index+1 == Gantt_chart_number[i][3]:
            
            
            each_Job  = str(Gantt_chart_number[i][0])
            
            st = Gantt_chart_number[i][1]  #改成Gantt_chart_number 開始時間
            # print(st)
            
            endt = Gantt_chart_number[i][2]  #改成Gantt_chart_number 結束時間
            # print(endt)
         
        
            df.append(dict(Task='%s'% ("Job"+ str(int(float(each_Job)))), Start='%s' % str(datetime.datetime.strptime("2021/11/14 00:00:00","%Y/%m/%d %H:%M:%S") + timedelta(seconds = st)),
                            Finish='%s'  % str(datetime.datetime.strptime("2021/11/14 00:00:00","%Y/%m/%d %H:%M:%S") +timedelta(seconds = endt)),Resource= value))
        
        
    return df

for index,value in enumerate(All_Machines_Name):
    # print(index,value)
    draw_machine(index ,value)
#------------------------------------

#呈現圖表
fig = px.timeline(df, x_start="Start", x_end="Finish",y="Resource",color="Task",text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
# fig = px.timeline(df, x_start="Start", x_end="Finish",y="Resource",color = "Task", hover_name = "Task",text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
fig.update_traces(textposition='inside',marker_line_color='rgb(8,48,107)')
# fig.update_layout(yaxis={'categoryorder':'category ascending'})  #將Y軸由小排到大
#將Y軸由大排到小
fig.update_layout(yaxis={'categoryorder':'category ascending'}, title={   # 设置整个标题的名称和位置
        "text":"甘特圖",
        "y":0.96,  # y轴数值
        "x":0.5,  # x轴数值
        "xanchor":"center",  # x、y轴相对位置
        "yanchor":"top"  
    })



fig.write_html("tsmc.html")#寫進網頁

plot(fig)
