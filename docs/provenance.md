# 来源与复用评估

检索日期：2026-09-28。实现依据为用户提供的 PDF：Süli & Trautwein, *Convergent numerical schemes for the viscoelastic Giesekus model in two dimensions*, arXiv:2512.22831v1, 2025-12-28。用户文件名的 2026 不改变 PDF 内版本。本次锁定 v1，不自动替换为在线新版本。

- 论文：https://arxiv.org/abs/2512.22831 。采用 (1.1)、(3.3)、(3.8)、(3.19) 和 5.1–5.2 节。
- deal.II 9.7.1：本地 `../Dealii/install-local`，作为链接依赖使用。许可证见其 `LICENSE.md`（主体 LGPL-2.1-or-later，部分文件另有双重许可）。未将 deal.II 源码复制到新仓库。
- 官方单纯形支持：https://dealii.org/9.6.0/doxygen/deal.II/group__simplex.html 。参考 FE_SimplexP、MappingFE、单纯形积分的 API 组合，并核对本地 9.7.1 头文件。
- 官方 step-57：https://dealii.org/9.2.0/doxygen/deal.II/step_57.html 。参考不可压 Navier–Stokes/Newton 的工程组织；其稳态/几何/离散细节不等同本论文，未直接复制为本算法。
- 官方开发仓库：https://github.com/dealii/dealii 。使用本地安装，不追踪 master。

检索了论文题目、作者及 Giesekus/deal.II/GitHub 组合词；未发现可确认与本论文 F 变量、时间层和稳定化一致的现成实现。相关 B 张量、薄膜近似及固体黏弹性实现不能直接替代 (3.3)。这一结论限于本次检索，不代表不存在作者私有或未检索到的代码。

新仓库中的算法代码自行编写，参考以上数学与 API 文档。尚未为新代码选择公开分发许可证，也未创建远程仓库或发布代码；未来公开时再由用户决定许可证。
