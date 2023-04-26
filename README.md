# single-cell
### 1. 常用分析软件
![Image](https://github.com/zhaohh52/single-cell/blob/main/support/softwares.jpg)
### 2. 如何理解 percent_mito,percent_ribo,percent_hb三个指标？   
      percent_hb(红细胞基因表达比例)：表明红细胞这个单细胞亚群的比例，一般来说不研究红细胞，所以过滤它没有问题。   
      percent_mito(线粒体基因表达比例)：表明细胞状态，值过高可能是濒临死亡的细胞，同样，不能一概而论，有些组织样本的细胞处于高代谢过程，该值会高于正常组织。   
      percent_ribo(核糖体基因表达比例）：我们之所以过滤这些，是因为在实际的实验操作过程中，会产生一些细胞杂质的影响，从而导致数据不准确。   
