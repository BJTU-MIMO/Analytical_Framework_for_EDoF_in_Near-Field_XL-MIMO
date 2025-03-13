# Analytical_Framework_for_Effective_Degrees_of_Freedom_in_Near-Field_XL-MIMO

This is a code package is related to the following scientific article:

Z. Wang, J. Zhang, W. Yi, H. Xiao, H. Du, D. Niyato, B. Ai, and D. W. K. Ng, "Analytical Framework for Effective Degrees of Freedom in Near-Field XL-MIMO," in IEEE Transactions on Wireless Communications, to appear, 2025.

Available at: https://arxiv.org/abs/2401.15280

The package contains a simulation environment, based on Matlab, that reproduces the numerical results in the article. *We encourage you to also perform reproducible research!*

## Abstract of Article
Extremely large-scale multiple-input-multiple-output (XL-MIMO) is an emerging transceiver technology for enabling next-generation communication systems, due to its potential for substantial enhancement in both the spectral efficiency and spatial resolution. However, the achievable performance limits of various promising XL-MIMO configurations have yet to be fully evaluated, compared, and discussed. In this paper, we develop an effective degrees of freedom (EDoF) performance analysis framework specifically tailored for near-field XL-MIMO systems. We explore five representative distinct XL-MIMO hardware designs, including uniform planar array (UPA)-based with infinitely thin dipoles, two-dimensional (2D) continuous aperture (CAP) plane-based, UPA-based with patch antennas, uniform linear array (ULA)-based, and one-dimensional (1D) CAP line segment-based XL-MIMO systems. Our analysis encompasses two near-field channel models: the scalar and dyadic Green's function-based channel models. More importantly, when applying the scalar Green's function-based channel, we derive EDoF expressions in the closed-form, characterizing the impacts of the physical size of the transceiver, the transmitting distance, and the carrier frequency. In our numerical results, we evaluate and compare the EDoF performance across all examined XL-MIMO designs, confirming the accuracy of our proposed closed-form expressions. Furthermore, we observe that with an increasing number of antennas, the EDoF performance for both UPA-based and ULA-based systems approaches that of 2D CAP plane and 1D CAP line segment-based systems, respectively. Moreover, we unveil that the EDoF performance for near-field XL-MIMO systems is predominantly determined by the array aperture size rather than the sheer number of antennas. 

## Content of Code Package

The package includes the source codes applied in this paper.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article:

```
@ARTICLE{10856805,
  author={Wang, Zhe and Zhang, Jiayi and Yi, Wenhui and Xiao, Huahua and Du, Hongyang and Niyato, Dusit and Ai, Bo and Ng, Derrick Wing Kwan},
  journal={IEEE Trans. Wireless Commun.}, 
  title={Analytical Framework for Effective Degrees of Freedom in Near-Field XL-MIMO}, 
  year={2025},
  month={to appear,},
  volume={},
  number={},
  pages={1-1}}
```
