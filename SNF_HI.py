...
...
#待论文发表后公开模型细节
"""
从这里开始进入了SNF相似性融合网络的步骤，上面使用MAD计算中位数偏差挑出了异质性相关的放射组学特征；
下面使用SNF聚类算法，得出所有病例的异质性聚类标签；
最后还需要根据MAD选出的这个放射组学特征的定义得出它们与异质性的正负相关关系，然后计算HI(异质性指数，越高代表异质性越高)；
HI的定义：经过Z-Score标准化处理之后的正相关特征与（负相关特征的相反数）的算术和。
"""
Sum_T1WDJY_T2W_features = []
Sum_T1WDJY_T2W_features.append(T1WD_features.values)
Sum_T1WDJY_T2W_features.append(T1WJ_features.values)
Sum_T1WDJY_T2W_features.append(T1WY_features.values)
Sum_T1WDJY_T2W_features.append(T2W_features.values)
Sum_T1WDJY_T2W_features.append(T1WD_features_dilation.values)
Sum_T1WDJY_T2W_features.append(T1WJ_features_dilation.values)
Sum_T1WDJY_T2W_features.append(T1WY_features_dilation.values)
Sum_T1WDJY_T2W_features.append(T2W_features_dilation.values)

# 输出相似性网络，可以通过调节mu的值来影响最终结果
# euclidean是欧几里得距离，sqeuclidean是平方欧几里得距离。K是邻近度的个数每个数据点的最近邻个数。mu是平滑参数，用于调整相似度矩阵的构建
# 设置随机种子
np.random.seed(2024)
affinity_networks = snf.make_affinity(Sum_T1WDJY_T2W_features, metric='sqeuclidean', K=13, mu=0.6)
# 融合网络,t是迭代次数
fused_network = snf.snf(affinity_networks,  t=20)
# 输出最优和次优的聚类族数目
best, second_best = snf.get_n_clusters(fused_network, n_clusters=range(2, 6))
print(best, second_best)
labels = spectral_clustering(fused_network, n_clusters=2)
print(labels)

"""
计算HI异质性分数:
步骤：1. 计算每个特征与标签的相关系数
      2. 分别输出正相关和负相关的特征
      3. 计算HI值：经过Z-Score标准化处理之后的正相关特征与（负相关特征的相反数）的算术和。
"""
#
Sum_T1WDJY_T2W_features_new =pd.concat([T1WD_features,T1WJ_features,T1WY_features,
                                    T2W_features,T1WD_features_dilation,T1WJ_features_dilation,
                                    T1WY_features_dilation,T2W_features_dilation],axis=1)  # 左右拼接
print(Sum_T1WDJY_T2W_features_new.shape)


Sum_T1WDJY_T2W_features_ZScore = StandardScaler().fit_transform(Sum_T1WDJY_T2W_features_new)
Sum_T1WDJY_T2W_features_ZScore_df = pd.DataFrame(Sum_T1WDJY_T2W_features_ZScore)
Sum_T1WDJY_T2W_features_ZScore_df.columns = (T1WD_features.columns.tolist() + T1WJ_features.columns.tolist() +
                                             T1WY_features.columns.tolist() + T2W_features.columns.tolist() +
                                             T1WD_features_dilation.columns.tolist() + T1WJ_features_dilation.columns.tolist() +
                                             T1WY_features_dilation.columns.tolist() + T2W_features_dilation.columns.tolist())

# 假设你的特征数据在df中，标签数据在labels中
df = Sum_T1WDJY_T2W_features_ZScore_df  # 示例数据
labels = pd.Series(labels)  # 示例标签
# # 计算每个特征与标签的相关系数
correlations = df.corrwith(labels)

# 分别输出正相关和负相关的特征名
positive_corr = correlations[correlations > 0].index.tolist()
negative_corr = correlations[correlations < 0].index.tolist()
# print("正相关的特征列索引:", positive_corr)
# print("负相关的特征列索引:", negative_corr)

Sum_T1WDJY_T2W_features_ZScore_df_positive = Sum_T1WDJY_T2W_features_ZScore_df[positive_corr]
Sum_T1WDJY_T2W_features_ZScore_df_negative = Sum_T1WDJY_T2W_features_ZScore_df[negative_corr]

# 计算HI值
HI = (Sum_T1WDJY_T2W_features_ZScore_df_positive.sum(axis=1) - Sum_T1WDJY_T2W_features_ZScore_df_negative.sum(axis=1))
HI_Zscored = (HI - HI.mean()) / HI.std()  # 标准化处理

Patient_Name_info = pd.read_excel("../Tumor_Heterogeneity_MRZ/Radiomics_feature_excels/Seclected_41_features/features_41_selected_and_patientsName.xlsx")
Patient_Name_info = Patient_Name_info['CureNo'].values.tolist()

HI_info = np.concatenate((np.array(Patient_Name_info).reshape(-1,1),np.array(labels).reshape(-1,1),
                          np.array(HI).reshape(-1,1), np.array(HI_Zscored).reshape(-1,1)), axis=1)  # 左右拼接
HI_df = pd.DataFrame(HI_info)
HI_df.columns = ['CureNo','SNF_Label', 'HI', 'HI_Zscored']
HI_df.to_excel("./MRZ_Excels/HI_df.xlsx",index=False)  # 输出HI值表格

...
...
