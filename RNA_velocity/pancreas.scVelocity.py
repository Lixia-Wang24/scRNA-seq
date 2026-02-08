#!/usr/bin/env python
# encoding: utf-8

#!/usr/bin/env python
# encoding: utf-8

import scvelo as scv
import os

# 把核心逻辑封装到一个函数里，或者直接放在 if 下面
def main():
    # 设置画图参数
    scv.settings.verbosity = 3
    scv.settings.presenter_view = True

    print(">>> 1. 正在加载 Pancreas 示例数据...")
    adata = scv.datasets.pancreas()

    print(">>> 2. 正在进行预处理 (Preprocessing)...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    print(">>> 3. 正在计算 RNA 速率 (Velocity)...")
    scv.tl.velocity(adata)

    # 【关键点】这一步使用了多进程，必须在 main 保护下运行
    print(">>> 3.1 计算速度图 (Velocity Graph)...")
    scv.tl.velocity_graph(adata)

    print(">>> 4. 正在保存结果...")
    adata.write('result_pancreas.h5ad')

    # 保存图片
    scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters', save='pancreas_embedding.png')

    print("✅ 分析完成！")
    print("数据已保存为: result_pancreas.h5ad")
    print("图片已保存到: figures/pancreas_embedding.png")

# 【这是您缺少的部分】
# 只有当这个文件被直接运行时，下面的代码才会被执行
if __name__ == "__main__":
    main()
