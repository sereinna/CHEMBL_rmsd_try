import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('results.csv', header=None)
data_array = data.dropna().values  # 删除所有包含NaN的行或列

# 筛选需要绘制的数据
mask = (data_array <= 600)  # 设置筛选条件
data_array = data_array[mask]

# 划分区间并计算每个区间内值的出现次数
num_bins = 2000
hist, bins = np.histogram(data_array, bins=num_bins)

# 计算每个区间的中点和宽度
bin_width = bins[1] - bins[0]
bin_centers = (bins[1:] + bins[:-1]) / 2

# 绘制直方图
plt.bar(bin_centers, hist, width=bin_width, alpha=0.5)

# 设置x轴标题和y轴标题
plt.xlabel("Value")
plt.ylabel("Frequency")

# 设置y轴范围
plt.ylim(0, np.max(hist) * 1.1)

plt.savefig("0_600_fenbu.svg", format="svg")

# 显示图像
plt.show()

--------------------------------------------------------------------------------------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('results.csv', header=None)
data_array = data.dropna().values  # 删除所有包含NaN的行或列
plt.rcParams.update({'font.size': 7})
# 设置子图的行数和列数
rows, cols = 3, 3

# 计算图像的总数
n = rows * cols

# 遍历每一列数据，并生成直方图
for i in range(n):
    if i >= data_array.shape[1]:
        break  # 如果列数超过数据的实际列数，则跳出循环
    
    # 筛选小于等于600的数据
    mask = (data_array[:, i] <= 200)
    filtered_data = data_array[mask, i]
    
    # 创建子图
    plt.subplot(rows, cols, i+1)
    
    # 对当前列进行划分区间并计算每个区间内值的出现次数
    hist, bins = np.histogram(filtered_data, bins=200)
    
    # 计算每个区间的中点和宽度
    bin_width = bins[1] - bins[0]
    bin_centers = (bins[1:] + bins[:-1]) / 2
    
    # 绘制直方图
    plt.bar(bin_centers, hist, width=bin_width, alpha=0.5)
    
    # 设置子图标题和x轴、y轴标题
    plt.title(f"Column {i+1}")
    plt.xlabel("Value")
    plt.ylabel("Frequency")

# 调整子图之间的距离和边距
plt.subplots_adjust(hspace=0.65, wspace=0.65)

# 保存图像为向量图
plt.savefig("histograms.svg", format="svg")

# 显示图像
plt.show()