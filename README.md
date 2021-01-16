# Ground-motion-generation-using-CWT
## 说明：

基于连续小波变换根据目标反应谱修正种子地震动得到目标地震动。

代码修改于：Montejo, L. A., & Suarez, L. E. (2013). An improved CWT-based algorithm for the generation of spectrum-compatible records. International Journal of Advanced Structural Engineering, 5(1), 26.  https://doi.org/10.1186/2008-6695-5-26 

## 使用方法：

1. 种子地震动文件：".//groundmotion/EW.txt"
2. 目标反应谱文件：".//TargetSpectrum//1//PredictedGM.txt"
3. 需要修改MainFuction.m中相应文件路径