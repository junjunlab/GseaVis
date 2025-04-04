
# GseaVis <img src="man/gseaVis-logo.png" align="right" height="200" />

<!-- badges: start -->

The goal of GseaVis is to visualize GSEA enrichment results as an implement package for **enrichplot** _gseaplot2_ function. And some codes origin from **enrichplot** package, thanks for **Guangchuang Yu** professor's contribution! The enrichment results from [**clusterProfiler**](https://github.com/YuLab-SMU/clusterProfiler) and [**GSEA software**](http://www.gsea-msigdb.org/gsea/index.jsp) can be supported as input for **GseaVis** for visualization.

You can mark your gene name on GSEA plot and this package also support more avaliable parameters to customize your own plot.

<!-- badges: end -->

## Installation

You can install the development version of GseaVis from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("junjunlab/GseaVis")
```

## Citation

> Jun Zhang, Hongyuan Li, Wenjun Tao, Jun Zhou. [GseaVis: An R Package for Enhanced Visualization of Gene Set Enrichment Analysis in Biomedicine](https://onlinelibrary.wiley.com/doi/full/10.1002/mdr2.70000). Med Research, 2025.

## Examples

![image](https://user-images.githubusercontent.com/64965509/198213474-43775942-1a40-4603-b1f8-2c2b0f2778e7.png)

## Documentation

**https://junjunlab.github.io/gseavis-manual/**

## Related blogs

> - [**GseaVis 优雅的可视化 GSEA 富集结果**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247501276&idx=1&sn=dce53570ae507affd283ade6bf13e635&chksm=c184ffadf6f376bb877733ac98f1bae3dbe3d1f1e019d9dc044e976dfc0f4197d6df832ea074&token=503374955&lang=zh_CN#rd)
> - [**GseaVis 的一些新功能**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247503821&idx=1&sn=452994f7744ef4ae9b0a84cfc82be016&chksm=c184f5bcf6f37caa2b30f5994e63ccf451f16e1f9e75db5131bc6004b3d0d91db048a44a0434&token=503374955&lang=zh_CN#rd)
> - [**GseaVis 对 gseKEGG 结果的支持优化**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247504498&idx=1&sn=9397b6e0ba0142e73df648bb86486003&chksm=c184e803f6f36115ccaf03dc792886a8cc5ea27e7b6c6e78282c022ad57bc73a7ee00d60d0e7&token=503374955&lang=zh_CN#rd)
> - [**R 包 bug 修复及问题一览**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247506202&idx=1&sn=e0e464ea398b5f53660109dc1bb82ea5&chksm=c184e36bf6f36a7d86509fd6a691c8e421b9143867664f0c18b5b010d558f4f37eccd2992aa9&token=503374955&lang=zh_CN#rd)
> - [**dotplotGsea 可视化 GSEA 富集结果**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247506271&idx=1&sn=74c8eabda17915d1931bc6fec7e054f2&chksm=c184e32ef6f36a38974b53c7586a00d192753e691be0ce09cfe501e4ea493ec1c2bbfbce25bc&token=503374955&lang=zh_CN#rd)
> - [**GSEA 结果火山图可视化**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247506419&idx=1&sn=f9ad91a426fc6ef1f6fe7e5ce0343845&chksm=c184e382f6f36a94c5ac6cdd662eb69ef6736d778c64635ada874d6f512f5cd91935047cfd72&token=503374955&lang=zh_CN#rd)
> - [**GseaVis 多条通路可视化**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247506468&idx=1&sn=784d933d674ebaa7eb47cc02597a8437&chksm=c184e055f6f3694356c618a07299d6215c1a87521bab4553b2b34938ae7799d1bd9927306fa2&token=503374955&lang=zh_CN#rd)
> - [**GseaVis 为你的通路基因添加热图**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247506582&idx=1&sn=24994c3eb73d5c30d5c56ae4894151a6&chksm=c184e0e7f6f369f151009e47e59a0ddbb7727dcd4dce8503f8fe861f1cf08103439b791a3eec&token=503374955&lang=zh_CN#rd)
> - [**bugs 报告和修复**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247507735&idx=1&sn=d8236c12a07beecc5d6c181b196a9a78&chksm=c184e566f6f36c7072f382be27259127b4fa9c0b1228c891f5cfc35869861b3d9b8f6e9b0824&token=139164705&lang=zh_CN#rd)
> - [**同学你又在画 GSEA?**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247507943&idx=1&sn=2dc950650892f93a53eb1ef9abab6555&chksm=c184e596f6f36c8073acb3b614897062249c9b51c3aff4657dc9156670a9ff7a99df78dc7787&token=495330596&lang=zh_CN#rd)
> - [**GseaVis 一键对接 GSEA 软件结果并可视化**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247508126&idx=1&sn=99fb4220166f8865762a6c2eb495ebe4&chksm=c1849aeff6f313f98b412186a5139e72ebb147f1899c00e7919899be6849f604667d987a1090&token=1432898004&lang=zh_CN#rd)
> - [**听说你想把多个样本的 GSEA 画在一起?**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247509949&idx=1&sn=712115d6aba5e45c19da2905308cfe3b&chksm=c1849dccf6f314da0d6af67f840410256c03d7bcefe5a576bc7babe46422bdf61c8cc44b81b2&token=353110902&lang=zh_CN#rd)
> - [**听说你只有富集表格还想画 GSEA 图?**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247512646&idx=1&sn=d9241b243771dcb355c74fd1a7f9431c&chksm=c1848837f6f301216f9d4c06659f01814b54994b552b5f15d9df9c5b227f417288a7cd4061e6&token=1853819831&lang=zh_CN#rd)
> - [**GSEA 还可以环形可视化？**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247514139&idx=1&sn=eae730e19b0ba2fa1f90dfb3ac10c582&chksm=c184826af6f30b7caa9d07c496d8237d44e01c9440bea1d0ac8b73846871a05e178392d0c2dc&token=1396583615&lang=zh_CN#rd)
> - [**听说你想把 fgsea 结果用 GseaVis 美化？**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247514164&idx=1&sn=0b7278817ff0a773cdfd89d3e373576f&chksm=c1848245f6f30b53b8ed7138363e72ec8a882d34f454a463a7ce79b840153c91712045ef1f97&token=1396583615&lang=zh_CN#rd)
> - [**桑基图加富集图一行代码出图?**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247513953&idx=1&sn=7c9eb4545986fa9187cdd66cae8b1456&chksm=c1848d10f6f30406d0bc3a86b3aa0150eb7191aea24bc4b9e366e0c04bc0be13a0a382253dfa&token=96515316&lang=zh_CN#rd)
> - [**GseaVis 绘制山脊图**](https://mp.weixin.qq.com/s/yOB9kFWtNYAqMLHoZ765eg?token=1296025231&lang=zh_CN)
