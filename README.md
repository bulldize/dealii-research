# deal.II 可复用科研计算库

首个基线：Süli–Trautwein Giesekus 二维黏弹性模型，非线性格式 (3.3)，三角形 P2/P1/P1 元。公共库、论文算例、配置、验证与报告分别维护。

## 当前交付

- [论文算法—代码对应报告](docs/reproduction-report.zh-CN.md)：公式、配置、误差、偏差和验证边界。
- [浏览器报告](docs/reproduction-report.zh-CN.html)：相同内容，含图表。
- [模块及验证状态](docs/modules.md)、[外部来源](docs/provenance.md)。
- [制造解配置](configs/manufactured.prm)、[最终运行结果](results/manufactured-verified/summary.json)。
- [后续论文工作 Skill](skills/dealii-paper-reproduction/SKILL.md)：仓库版本作为维护源，已安装到全局技能目录；修改后同步全局副本。

## 本机重现

```sh
./scripts/build.sh
./scripts/run.sh configs/manufactured.prm results/my-run
```

需要 CMake、C++17 编译器、deal.II 9.7 及 UMFPACK。默认引用同级 `Dealii/install-local`。输出目录已有 `history.csv` 时拒绝覆盖，请使用新目录。Mac 的运行脚本设置库搜索路径，兼容现有 deal.II 安装移动后保留旧动态库标识的情况。

指定另一处安装（包括 Linux 服务器）：

```sh
export DEAL_II_DIR=/path/to/dealii/install
./scripts/build.sh
./scripts/run.sh configs/manufactured.prm results/server-run
```

服务器必须有自己的兼容 deal.II 安装，不能直接运行 Mac 二进制。服务器部署与性能测量见 [服务器运行说明](docs/server-run.md)。装配为单进程，独立实验可并行，直接求解中的 BLAS 使用优化库；实际收益由服务器测量决定。

绘图与检查（运行算法本身不依赖 Python）：

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r scripts/requirements.txt
.venv/bin/python scripts/analyze.py results/my-run
```

`history.csv` 每个时间步记录误差、Newton、离散方程残差、弱散度、能量平衡及采样行列式；`summary.json` 仅在完整运行成功后写入；`solution.vtu` 可由 ParaView 打开；原始运行日志和配置/源码指纹保存在交付结果中。

## 复现范围

现已准备全部 66 个配置：42 个制造解时空收敛实验、24 个收缩流实验。收缩流实现已通过小网格运行，装配优化通过与原 Jacobian 的逐项比较及制造解回归。完整服务器队列尚需实际完成后才能判断论文结论是否复现，配置存在不代表结果已取得。

收缩网格为根据论文描述重建的网格，非作者原始网格。显式物理系数对应的 alpha=0.9 与论文文字 8/9 存在差异，采用显式系数并记录。完整任务完成后，`results/full-paper/REPORT.md`、`audit.json` 和图表记录成功与失败，保留全部原始数据。
