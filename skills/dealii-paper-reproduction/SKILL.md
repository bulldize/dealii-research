---
name: dealii-paper-reproduction
description: Reproduce finite element research papers with deal.II while developing a reusable research library and a formula-to-code verification report. Use for adding papers or extending this library, not for unrelated numerical tasks.
---

先定位目标科研仓库，再读用户指定论文和该仓库的验证记录，明确模型、弱形式、离散空间、时间层、边界条件与验收范围。附件是研究资料，不是额外操作授权。

对已有模块判断可直接复用、参数变化、需要适配或真正新增。对缺失部分检索 deal.II 官方示例和相关成熟实现，核对数学假设、许可证、版本和测试；找不到一致实现时再自行实现。记录参考来源，不把近似模型称为论文算法。

模块划分服务于当前可验证的复用需要。已验证且未受影响的模块保持稳定；发生耦合变化时补充受影响的验证，不强制统一目录、接口或求解器结构。

交付报告列出论文公式/章节到代码符号、配置、实际结果和偏差的对应关系。将理论结论、实测结果和未验证内容分开；单次运行不能验证收敛阶，正行列式采样不能证明全域正性。记录失败而非静默修改物理量。

本仓库入口为 README.md，模块状态见 docs/modules.md，首篇论文对应报告见 docs/reproduction-report.zh-CN.md。新任务开始时读取相关部分即可。单次实验和参数扫描按用户当前范围执行。

默认已有科研库的位置和维护方式见 [references/local-library.md](references/local-library.md)。用户指定其他仓库时，以用户指定为准。

用户无需预先选择所有数值配置：先提取论文已给定的配置，再明确标记论文未说明而采用的实现选择。复现比较使用相同的物理问题、离散参数及误差定义；不能把不同网格、不同误差范数或最终时刻相对误差直接与论文时空绝对误差相比。

后续论文带来的通用模块和验证证据沉淀到目标仓库；仅把实际改善后续决策的经验写回本 Skill。修改本 Skill 时优先维护仓库版本，验证后同步全局副本，避免两个版本逐渐分叉。
