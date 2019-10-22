# CF-tracking-for-remote-sensing-image
## Introduction
一路靠着运气混进了决赛，最后混了个三等奖。实际上我们当时的思路也很简单，深度学习做跟踪这一块儿我们本身就不是很熟，所以想着用最简单的模型来做，所以就用了最简单的Standard DCF[1]。针对跟踪进度不高，我们做了三点改进，当然这三点改进都很简单，主要也是借鉴了别的算法的思想。  
1、训练阶段，增加一些样本，希望提高模型的判别能力，这主要是CA-CF[2]的思想，我们直接拿来用了。  
2、针对完全遮挡问题，用检测响应图来判断目标是否丢失，丢失后重新搜索。这实际上是长期跟踪LCT[3]的思想，我们图快做了简化，不用单独的分类器和滤波器，直接用跟踪的检测响应图，当然问题也很明显，就是参数必须手动给。  
3、针对目标的旋转，感觉这是最难的一部分。标准的相关滤波算法是没办法处理这个问题的（有一篇文章[4]倒是说了这个问题,不过没太看懂，也就不知道咋用），这种旋转带来的不仅是变化的尺寸长宽比，还有特征的变化。我们使用的是传统的手工特征，处理这种大范围的旋转跟踪经常丢失，像ECO那样用深度特征应该会好一点。最后的做法其实还是模板匹配的思想，不过改成了用稀疏表示分类[5]来做。将尺度搜索改为角度搜索，通过角度变化来估计尺度变化（其实算是投机取巧了，不过稀疏确实很好用啊）。  
## 参考文献和代码
感谢这些大佬，让我等门外汉也能一窥跟踪算法之妙  
[1] KCF/DCF: João F. Henriques, Rui Caseiro, Pedro Martins, Jorge Batista. "High-Speed Tracking with Kernelized Correlation Filters." TPAMI (2015).  
[2] CF+CA: Matthias Mueller, Neil Smith, Bernard Ghanem. "Context-Aware Correlation Filter Tracking." CVPR (2017).   
[3] LCT: Chao Ma, Xiaokang Yang, Chongyang Zhang, Ming-Hsuan Yang. "Long-term Correlation Tracking." CVPR (2015).   
[4] IBCCF: Feng Li, Yingjie Yao, Peihua Li, David Zhang, Wangmeng Zuo, Ming-Hsuan Yang. "Integrating Boundary and Center Correlation Filters for Visual Tracking With Aspect Ratio Variation." ICCV workshop (2017).   
[5] L1APG: C. Bao, Y. Wu, H. Ling and H. Ji, "Real time robust L1 tracker using accelerated proximal gradient approach", IEEE Conf. on Computer Vision and Pattern Recognition (CVPR), Rhode Island, 2012.   
